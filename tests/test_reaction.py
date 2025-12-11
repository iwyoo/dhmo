import os
import pytest
from src.dhmo import load


def test_dehydration_reaction():
    # Example yaml path
    yaml_path = os.path.join(os.path.dirname(__file__), '../examples/condensation.yaml')
    reaction = load(yaml_path)
    rxn_smarts = reaction.rxn_smarts
    assert isinstance(rxn_smarts, str)
    assert '>>' in rxn_smarts
    print(f"Generated reaction SMARTS: {rxn_smarts}")

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