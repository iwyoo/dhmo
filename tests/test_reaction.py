import os
import pytest
from src.dhmo import load


def test_dehydration_reaction():
    # Example yaml path
    yaml_path = os.path.join(os.path.dirname(__file__), '../examples/condensation.yaml')
    reaction = load(yaml_path)
    smarts = reaction.smarts
    assert isinstance(smarts, str)
    assert '>>' in smarts
    print(f"Generated SMARTS: {smarts}")

    input_smiles = ["CCO", "CC(=O)O"]
    output_smiles = reaction.run(input_smiles)
    print(f"Input SMILES list: {input_smiles}")
    print(f"Output SMILES list: {output_smiles}")