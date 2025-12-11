from itertools import product
from typing import Any, Dict

from rdkit import Chem
from rdkit.Chem import AllChem


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
        self.rxn_smarts = self._preprocess_smarts(
            self.input_smarts_map,
            self.output_smarts_map,
        )
        self.rxn = AllChem.ReactionFromSmarts(self.rxn_smarts)

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

            smarts = self._apply_atom_mapping_by_symbol(pattern, mapping)
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

    def _apply_atom_mapping_by_symbol(self, pattern: str, mapping: str) -> str:
        """
        Convert a SMILES pattern and atom mapping to a SMARTS string with atom mapping.
        For example, given the SMILES pattern "C(=O)O" and mapping "1__2__",
        the function returns "[C:1](=[O:2])O".

        Args:
            pattern (str): SMARTS pattern (e.g., "C(=O)O")
            mapping (str): Atom mapping string (e.g., "1__2__")
        Returns:
            str: SMARTS string with atom mapping (e.g., "[C:1](=[O:2])O")
        """
        mol = Chem.MolFromSmarts(pattern)
        if mol is None:
            raise ValueError(f"Invalid SMILES pattern: {pattern}")

        pattern_lower = pattern.lower()
        atom_pos2index = {}
        i = 0
        for atom_idx, atom in enumerate(mol.GetAtoms()):
            symbol = atom.GetSymbol()
            i = pattern_lower.find(symbol.lower(), i)
            if i == -1 or i >= len(mapping):
                raise ValueError(
                    f"Mapping does not match pattern for atom '{symbol}' in pattern '{pattern}' with mapping '{mapping}'"
                )

            if len(symbol) == 2 and pattern[index + 1] != "_":
                raise ValueError(
                    f"Second character of atom symbol '{symbol}' in pattern '{pattern}' must be followed by an underscore in mapping"
                )
            atom_pos2index[i] = atom_idx
            i += 1

        map_num2pos = {num: pos for pos, num in enumerate(mapping) if num != "_"}

        if len(map_num2pos) > len(atom_pos2index):
            raise ValueError(
                f"More mapped atoms in mapping '{mapping}' than atoms in pattern '{pattern}'"
            )

        for num, pos in map_num2pos.items():
            atom_idx = atom_pos2index.get(pos)
            if atom_idx is None:
                raise ValueError(
                    f"No atom found at position {pos} in pattern '{pattern}' for mapping '{mapping}'"
                )
            atom = mol.GetAtomWithIdx(atom_idx)
            atom.SetAtomMapNum(int(num))

        smarts = Chem.MolToSmarts(mol)
        return smarts

    def _generate_reactant_combinations(self, input_mols: list) -> list:
        """
        Generate all possible combinations of reactant molecules for the reaction.
        Args:
            input_mols (list of rdkit.Chem.Mol): List of input molecule objects.
        Returns:
            list of tuple: List of tuples, each containing a combination of reactant molecules.
        """
        possible_mol_lists = []
        for input_smarts in self.input_smarts_map.values():
            possible_input_mols = []
            pattern_mol = Chem.MolFromSmarts(input_smarts)
            for mol in input_mols:
                if mol.HasSubstructMatch(pattern_mol):
                    possible_input_mols.append(mol)
            possible_mol_lists.append(possible_input_mols)
        combinations = list(product(*possible_mol_lists))
        return combinations

    def run(self, inputs: list[str], ordered: bool = False):
        """
        Run the reaction on the given input molecules.
        Args:
            inputs (list of str): Input molecule SMILES list.
            ordered (bool): If True, use input order strictly. If False, try all combinations.
        Returns:
            list of str: Output molecule SMILES list.
        """
        if not ordered:
            # Remove duplicate inputs for unordered mode
            inputs = list(set(inputs))

        try:
            input_mols = [Chem.MolFromSmiles(smi) for smi in inputs]
        except Exception as e:
            raise ValueError(f"Error parsing input SMILES: {e}")

        if any(mol is None for mol in input_mols):
            raise ValueError("One or more input SMILES are invalid.")

        outputs = []
        if ordered:
            # Only one combination: input order as given
            combinations = [tuple(input_mols)]
        else:
            combinations = self._generate_reactant_combinations(input_mols)

        for mol_combination in combinations:
            self.rxn.Initialize()
            products_sets = self.rxn.RunReactants(mol_combination)
            for products in products_sets:
                output_combination = []
                is_product_valid = True
                for product in products:
                    try:
                        Chem.SanitizeMol(product)
                    except Exception:
                        is_product_valid = False
                        break
                    smi = Chem.MolToSmiles(product, isomericSmiles=True)
                    output_combination.append(smi)
                if is_product_valid:
                    outputs.append(".".join(output_combination))
        outputs = list(set(outputs))
        return outputs