import yaml

from .reaction import Reaction


_NAME_KEY = "name"
_DESCRIPTION_KEY = "description"
_INPUT_KEY = "inputs"
_OUTPUT_KEY = "outputs"
_CONDITION_KEY = "conditions"


def load(path):
    with open(path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)
    return Reaction(**data)