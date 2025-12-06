import yaml
import os
import pytest
from src.dhmo.reaction import Reaction


def load_yaml(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)


def test_dehydration_reaction():
    data = load_yaml(os.path.join(os.path.dirname(__file__), '../examples/condensation.yaml'))
    reaction_info = data['reaction']
    reaction = Reaction(
        name=reaction_info.get('name', ''),
        description=reaction_info.get('description', ''),
        inputs=reaction_info['inputs'],
        outputs=reaction_info['outputs'],
    )
    smarts = reaction.smarts
    assert isinstance(smarts, str)
    assert '>>' in smarts
    print(f"Generated SMARTS: {smarts}")