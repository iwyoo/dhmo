# DHMO (*/dɯ.mo/*)

DHMO is a Python library inspired by the famous 'dihydrogen monoxide' chemistry joke. DHMO focuses on making chemical reaction modeling as human-readable and accessible as possible, providing a friendly wrapper for defining and simulating chemical reactions.

## Features
- Human-friendly YAML-based chemical reaction definitions
- Easy-to-use reaction objects for simulating chemical transformations
- Focus on clarity and accessibility for chemists and non-experts alike

## Example: Condensation Reaction (YAML)

Below is an example YAML file describing a condensation reaction (see `examples/condensation.yaml`):

```yaml
name: Condensation
description: "Condensation reaction between an alcohol and a carboxylic acid to form an ester and water"
inputs:
  alcohol:
    - "*O"
    - "1_"
  acetic acid:
    - "*(=O)O"
    - "2_____"
outputs:
  ester:
    - "*(=O)O*"
    - "1_____2"
  water:
    - "O"
    - "_"
```

## Usage Example

You can load a chemical reaction from a YAML file and run it on input molecules as follows:

```python
import dhmo

# Load the reaction from YAML
yaml_path = "examples/condensation.yaml"
reaction = dhmo.load(yaml_path)

# Run the reaction with input SMILES
input_smiles_list = ["CCO", "CC(=O)O"]
output_smiles_list = reaction.run(input_smiles_list)
print(f"Input SMILES: {input_smiles_list}")
print(f"Output SMILES: {output_smiles_list}")
```
