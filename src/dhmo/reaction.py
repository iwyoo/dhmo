from rdkit import Chem
from rdkit.Chem import AllChem

from typing import Any, Dict


class Reaction:
    def __init__(
        self,
        *args: Any,
        name: str = "",
        description: str = "",
        inputs: Dict[str, list],
        outputs: Dict[str, list],
        **kwargs: Any
    ) -> None:
        self.name = name
        self.description = description
        self.input_smarts_map = self._preprocess_reaction_members(inputs)
        self.output_smarts_map = self._preprocess_reaction_members(outputs)
        self.smarts = self._preprocess_smarts(
            self.input_smarts_map,
            self.output_smarts_map,
        )

    def _preprocess_reaction_members(self, members: Dict[str, list]) -> Dict[str, str]:
        preprocessed_members = {}
        for member, rule in members.items():
            if not isinstance(member, str):
                raise ValueError(f"Member key must be a string, got: {repr(member)}")

            if not isinstance(rule, list):
                raise ValueError(
                    f"Rule for member '{member}' must be a list, got: {repr(rule)}"
                )

            if len(rule) != 2:
                raise ValueError(
                    f"Rule for member '{member}' must have exactly 2 elements, got: {repr(rule)}"
                )

            if not all(isinstance(rpart, str) for rpart in rule):
                raise ValueError(
                    f"All elements in rule for member '{member}' must be strings, got: {repr(rule)}"
                )

            pattern, mapping = rule
            if len(pattern) != len(mapping):
                raise ValueError(
                    f"Pattern and mapping for member '{member}' must have the same length, got: pattern={repr(pattern)}, mapping={repr(mapping)}"
                )
            
            if not all((c == "_" or c.isdigit() and c != '0') for c in mapping):
                raise ValueError(
                    f"Mapping for member '{member}' must contain only digits or underscores, got: {repr(mapping)}"
                )

            smarts = self._convert_smiles_to_smarts(pattern, mapping)
            preprocessed_members[member] = smarts
        return preprocessed_members

    @staticmethod
    def _preprocess_smarts(input_members: Dict[str, str], output_members: Dict[str, str]) -> str:
        """
        Preprocess the reaction SMARTS string from input and output member SMARTS maps.

        Args:
            input_members (dict): Mapping of input member names to SMARTS strings.
            output_members (dict): Mapping of output member names to SMARTS strings.
        Returns:
            str: Reaction SMARTS string.
        """
        input_smarts_part = []
        for input_member in input_members.values():
            input_smarts_part.append(f"({input_member})")
        input_smarts = ".".join(input_smarts_part)

        output_smarts_part = []
        for output_member in output_members.values():
            output_smarts_part.append(f"({output_member})")
        output_smarts = ".".join(output_smarts_part)

        reaction_smarts = f"{input_smarts}>>{output_smarts}"
        return reaction_smarts

    def _convert_smiles_to_smarts(self, pattern: str, mapping: str) -> str:
        """
        Convert a SMILES pattern and atom mapping to a SMARTS string with atom mapping.
        For example, given the SMILES pattern "C(=O)O" and mapping "1__2__",
        the function returns "[C:1](=[O:2])O". 

        Args:
            pattern (str): SMILES pattern (e.g., "C(=O)O")
            mapping (str): Atom mapping string (e.g., "1__2__")
        Returns:
            str: SMARTS string with atom mapping (e.g., "[C:1](=[O:2])O")
        """
        mol = Chem.MolFromSmiles(pattern)
        if mol is None:
            raise ValueError(f"Invalid SMILES pattern: {pattern}")
        
        # Assign each atom its mapping number
        pattern_lower = pattern.lower()
        i = 0
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()

            # Find the corresponding position in the pattern
            index = pattern.find(symbol, i)
            if index == -1 or index >= len(mapping):
                raise ValueError(
                    f"Mapping does not match pattern for atom '{symbol}' in pattern '{pattern}' with mapping '{mapping}'"
                )
            
            if len(symbol) == 2 and pattern[index + 1] != "_":
                raise ValueError(
                    f"Second character of atom symbol '{symbol}' in pattern '{pattern}' must be followed by an underscore in mapping"
                )

            if mapping[index] == "_":
                atom.SetAtomMapNum(0)
            else:
                atom.SetAtomMapNum(int(mapping[index]))
            i += 1

        smarts = Chem.MolToSmarts(mol)
        return smarts

    def run(self, inputs: list[str]):
        """
        Run the reaction on the given input molecules.

        Args:
            inputs (list of str): Input molecule SMILES list.
        Returns:
            list of str: Output molecule SMILES list.
        """

        rxn = AllChem.ReactionFromSmarts(self.smarts)
        try:
            input_mols = [Chem.MolFromSmiles(smi) for smi in inputs]
        except Exception as e:
            raise ValueError(f"Error parsing input SMILES: {e}")

        if any(mol is None for mol in input_mols):
            raise ValueError("One or more input SMILES are invalid.")

        products_sets = rxn.RunReactants(input_mols)
        outputs = []
        for products in products_sets:
            output_combination = []
            for product in products:
                smi = Chem.MolToSmiles(product, isomericSmiles=True)
                output_combination.append(smi)
            outputs.append(".".join(output_combination))
        return outputs