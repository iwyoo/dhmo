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

    input_smiles = ["CCO", "CC(=O)O"]
    output_smiles = reaction.run(input_smiles)
    print(f"Input SMILES list: {input_smiles}")
    print(f"Output SMILES list: {output_smiles}")

    # Reverse reaction test
    input_smiles = ["CC(=O)O", "CCO"]
    output_smiles = reaction.run(input_smiles)
    print(f"Input SMILES list: {input_smiles}")
    print(f"Output SMILES list: {output_smiles}")

    output_smiles = reaction.run(input_smiles, ordered=True)
    if len(output_smiles) > 0:
        raise AssertionError("Expected no products for ordered reaction with given inputs.")