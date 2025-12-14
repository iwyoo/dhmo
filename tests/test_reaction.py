import os
import pytest
from src.dhmo import load


def test_condensation_reaction():
    # Example yaml path
    yaml_path = os.path.join(os.path.dirname(__file__), "../examples/condensation.yaml")
    reaction = load(yaml_path)

    # Test running the reaction (unordered)
    input_smiles = ["CCO", "CC(=O)O"]
    output_smiles = reaction.run(input_smiles)
    if output_smiles != ["CC(=O)OC(C)=O.O", "CCOC(C)=O.O"]:
        raise AssertionError(f"Unexpected products: {output_smiles}")

    output_smiles = reaction.run(input_smiles, ordered=True)
    if output_smiles != ["CCOC(C)=O.O"]:
        raise AssertionError(f"Unexpected ordered products: {output_smiles}")

    # Test reversed reaction
    reversed_reaction = reaction.reversed()
    output_smiles_list = reversed_reaction.run(["CCOC(C)=O.O"])
    if output_smiles_list != ["CC(=O)O.CCO"]:
        raise AssertionError(f"Unexpected products from reversed reaction: {output_smiles_list}")

    # Test atom mapping preservation
    input_smiles = ["[CH3:1][CH2:2][OH:3]", "[CH3:4][C:5](=[O:6])[OH:7]"]
    output_smiles = reaction.run(input_smiles, ordered=True)
    if output_smiles != ["[CH3:1][CH2:2][O:3][C:5]([CH3:4])=[O:6].[OH2:7]"]:
        raise AssertionError(f"Unexpected products with atom-mapping: {output_smiles}")


def test_oxidation_of_primary_alcohol():
    yaml_path = os.path.join(os.path.dirname(__file__), "../examples/oxidation.yaml")
    reaction = load(yaml_path)

    # Ethanol oxidation: CCO (Expect acetaldehyde: CC=O)
    input_smiles = ["CCO"]
    output_smiles = reaction.run(input_smiles)

    if output_smiles != ["CC=O"]:
        raise AssertionError(f"Unexpected products: {output_smiles}")

    # Test atom index preservation for c1ccc([CH2:1][OH:2])cc1
    input_smiles = ["c1ccc([CH2:1][OH:2])cc1"]
    output_smiles = reaction.run(input_smiles)
    if output_smiles != ["c1ccc([CH:1]=[O:2])cc1"]:
        raise AssertionError(f"Unexpected products with atom-mapping: {output_smiles}")
