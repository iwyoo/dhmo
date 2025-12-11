import os
import pytest
from src.dhmo import load


def test_condensation_reaction():
    # Example yaml path
    yaml_path = os.path.join(os.path.dirname(__file__), '../examples/condensation.yaml')
    reaction = load(yaml_path)

    # Test running the reaction (unordered)
    input_smiles = ["CC(=O)O", "CCO"]
    output_smiles = reaction.run(input_smiles)
    if output_smiles != ['CCOC(C)=O.O']:
        raise AssertionError(f"Unexpected products: {output_smiles}")

    output_smiles = reaction.run(input_smiles, ordered=True)
    if len(output_smiles) > 0:
        raise AssertionError("Expected no products for ordered reaction with given inputs.")

    # Test reversed reaction
    reversed_reaction = reaction.reversed()
    output_smiles_list = reversed_reaction.run(["CCOC(C)=O.O"])
    if output_smiles_list != ['CC(=O)O.CCO']:
        raise AssertionError(f"Unexpected products from reversed reaction: {output_smiles_list}")


def test_oxidation_of_primary_alcohol():
    yaml_path = os.path.join(os.path.dirname(__file__), '../examples/oxidation.yaml')
    reaction = load(yaml_path)

    # Ethanol oxidation: CCO + [O] -> acetaldehyde + water
    input_smiles = ["CCO"]
    output_smiles = reaction.run(input_smiles)
    # Expect acetaldehyde (CC=O) and water (O)
    if output_smiles != ["CC=O"]:
        raise AssertionError(f"Unexpected products: {output_smiles}")
