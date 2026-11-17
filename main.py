
from pathlib import Path
from ase import Atoms
import numpy as np

# Extract data from the data file

def _parse_n2p2_block(block):
    symbols = []
    positions = []
    cell = []
    energy = None
    has_cell = False

    for line in block:
        parts = line.split()

        if not parts:
            continue

        if parts[0] == "lattice":
            cell.append([float(x) for x in parts[1:4]])
            has_cell = True

        elif parts[0] == "atom":
            positions.append([float(x) for x in parts[1:4]])
            symbols.append(parts[4])

        elif parts[0] == "energy":
            energy = float(parts[1])

    if energy is None:
        raise ValueError("Energy missing in n2p2 block")

    if not symbols:
        raise ValueError("No atoms found in n2p2 block")

    atoms = Atoms(
        symbols=symbols,
        positions=positions,
        cell=cell if has_cell else None,
        pbc=has_cell
    )

    return atoms, energy


def read_n2p2(path):
    atoms_list = []
    energies = []
    block = []
    inside_block = False

    with open(path, "r") as file:
        for line in file:
            line = line.strip()

            if line == "begin":
                block = []
                inside_block = True

            elif line == "end":
                if not inside_block:
                    raise ValueError("Unexpected end")

                atoms, energy = _parse_n2p2_block(block)
                atoms_list.append(atoms)
                energies.append(energy)
                inside_block = False

            elif inside_block:
                block.append(line)

    if inside_block:
        raise ValueError("Unclosed n2p2 block")

    return atoms_list, np.array(energies)


# Example usage
if __name__ == "__main__":
    path = Path(__file__).resolve().parent / "input.data"

    atoms_list, energies = read_n2p2(path)

    print("Number of structures:", len(atoms_list))
    print("Reference energies:", energies)
