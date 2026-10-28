import math

pdb_file = "data/raw/1BC8.pdb"

zn_coord = None
protein_atoms = []

with open(pdb_file, "r") as file:
    for line in file:

        if line.startswith("HETATM") and " ZN " in line:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            zn_coord = (x, y, z)

        elif line.startswith("ATOM"):
            chain = line[21]
            residue_name = line[17:20].strip()
            residue_number = line[22:26].strip()
            atom_name = line[12:16].strip()

            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            protein_atoms.append(
                (
                    chain,
                    residue_name,
                    residue_number,
                    atom_name,
                    x,
                    y,
                    z
                )
            )

print("Zn coordinate:", zn_coord)

cutoff = 5.0

nearby_residues = set()

for atom in protein_atoms:

    chain, residue_name, residue_number, atom_name, x, y, z = atom

    distance = math.sqrt(
        (x - zn_coord[0]) ** 2
        + (y - zn_coord[1]) ** 2
        + (z - zn_coord[2]) ** 2
    )

    if distance <= cutoff:
        nearby_residues.add(
            (chain, residue_name, residue_number)
        )

print("\nResidues within 5 Å of Zn:")

for residue in sorted(nearby_residues):
    print(residue)