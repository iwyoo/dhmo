import re
from typing import Dict, List

def parse_smarts_atoms(smarts: str) -> Dict[int, List[int]]:
    """
    Parse a SMARTS string and return a dictionary mapping atom indices to their positions in the string.
    Each key is the atom index (0-based), and the value is a list of positions
    (all character indices in the SMARTS string for that atom).
    Args:
        smarts (str): SMARTS string.
    Returns:
        Dict[int, List[int]]: {atom_idx: [positions in smarts string]}
    """
    # Atom in square brackets: \[[^\]]+\]
    # Atom outside square brackets: *, Br, Cl, N, O, S, P, F, I, b, c, n, o, s, p
    atom_pattern = re.compile(r'(\[[^\]]+\]|\*|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p)')
    positions = {}
    idx = 0
    for match in atom_pattern.finditer(smarts):
        start = match.start()
        end = match.end()
        # All character indices occupied by this atom
        positions[idx] = list(range(start, end))
        idx += 1
    return positions
