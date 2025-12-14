from itertools import product
from typing import Any, Dict
import copy

from rdkit import Chem
from rdkit.Chem import AllChem

from .utils.parser import parse_smarts_atoms


class Reaction:
    def __init__(
        self,
        *args: Any,
        name: str = "",
        description: str = "",
        inputs: Dict[str, list],
        outputs: Dict[str, list],
        is_reversed: bool = False,
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
        self.is_reversed = is_reversed

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

            smarts = self._apply_atom_mapping(pattern, mapping)
            preprocessed_members[member] = smarts
        return preprocessed_members

    def _preprocess_smarts(self, input_members: Dict[str, str], output_members: Dict[str, str]) -> str:
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

    def _apply_atom_mapping(self, pattern: str, mapping: str) -> str:
        """
        Apply atom mapping to a SMARTS string using a mapping string.
        Uses parse_smarts_atoms to get atom index to SMARTS positions.
        Args:
            pattern (str): SMARTS pattern (e.g., "C(=O)O")
            mapping (str): Atom mapping string (e.g., "1__2__")
        Returns:
            str: SMARTS string with atom mapping
        """
        mol = Chem.MolFromSmarts(pattern)
        if mol is None:
            raise ValueError(f"Invalid SMARTS pattern: {pattern}")

        atom_pos_dict = parse_smarts_atoms(pattern)

        # Collect all used positions from atom_pos_dict
        used_positions = set()
        for pos_list in atom_pos_dict.values():
            used_positions.update(pos_list)

        # Check unused positions in mapping before atom mapping
        for i, c in enumerate(mapping):
            if i not in used_positions and c != '_':
                raise ValueError(
                    f"Mapping character at position {i} ('{c}') "
                    "is not used for any atom and must be '_'."
                )

        # Now apply mapping
        for atom_idx, positions in atom_pos_dict.items():
            digit_chars = [mapping[pos] for pos in positions if mapping[pos].isdigit()]
            atom = mol.GetAtomWithIdx(atom_idx)
            if "0" in digit_chars:
                raise ValueError(f"Mapping digit '0' is not allowed for atom index {atom_idx} in mapping string: {digit_chars}")
            if len(digit_chars) > 1:
                raise ValueError(f"Multiple mapping digits found for atom index {atom_idx} in mapping string: {digit_chars}")
            elif len(digit_chars) == 1:
                atom.SetAtomMapNum(int(digit_chars[0]))
            else:
                atom.SetAtomMapNum(0)
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

        # Helper function to extract reactant map info for atom-mapping restoration
        def get_reactants_for_maps(rxn):
            res = {}
            for i in range(rxn.GetNumReactantTemplates()):
                r = rxn.GetReactantTemplate(i)
                for at in r.GetAtoms():
                    if at.GetAtomMapNum():
                        res[at.GetAtomMapNum()] = i
            return res

        outputs = []
        if ordered:
            # Only one combination: input order as given
            combinations = [tuple(input_mols)]
        else:
            combinations = self._generate_reactant_combinations(input_mols)

        for mol_combination in combinations:
            self.rxn.Initialize()
            products_sets = self.rxn.RunReactants(mol_combination)
            rnos = get_reactants_for_maps(self.rxn)
            for products in products_sets:
                output_combination = []
                is_product_valid = True
                for product in products:
                    try:
                        Chem.SanitizeMol(product)
                    except Exception:
                        is_product_valid = False
                        break

                    # Assign atom-mapping based on reactants
                    # using old_mapno to reactant index mapping
                    for product_atom in product.GetAtoms():
                        if product_atom.HasProp("old_mapno"):
                            mapno = int(product_atom.GetProp("old_mapno"))
                        else:
                            continue

                        reactant_mol = mol_combination[rnos[mapno]]
                        reactant_atom = reactant_mol.GetAtomWithIdx(
                            product_atom.GetIntProp("react_atom_idx"),
                        )
                        reactant_atom_map_num = reactant_atom.GetAtomMapNum()
                        product_atom.SetAtomMapNum(reactant_atom_map_num)

                    smi = Chem.MolToSmiles(product, isomericSmiles=True)
                    output_combination.append(smi)
                if is_product_valid:
                    output_smiles = ".".join(output_combination)
                    output_smiles = Chem.CanonSmiles(output_smiles)
                    outputs.append(output_smiles)
        outputs = sorted(set(outputs))
        return outputs

    def reversed(self, name=None, description=None):
        """
        Create a new Reaction instance with inputs and outputs swapped (retrosynthesis).
        Args:
            name (str, optional): Name for the reversed reaction. If None, no name is set.
            description (str, optional): Description for the reversed reaction. If None, no description is set.
        Returns:
            Reaction: New Reaction instance with reversed inputs/outputs.
        """
        if name is None:
            name = f"Reversed reaction of: {self.name}"
        if description is None:
            description = f"Reversed reaction of: {self.description}"

        new_rxn = Reaction.__new__(Reaction)

        new_rxn.name = name
        new_rxn.description = description
        new_rxn.input_smarts_map = copy.deepcopy(self.output_smarts_map)
        new_rxn.output_smarts_map = copy.deepcopy(self.input_smarts_map)
        new_rxn.rxn_smarts = new_rxn._preprocess_smarts(
            new_rxn.input_smarts_map,
            new_rxn.output_smarts_map,
        )
        new_rxn.rxn = AllChem.ReactionFromSmarts(new_rxn.rxn_smarts)
        new_rxn.is_reversed = not self.is_reversed
        return new_rxn