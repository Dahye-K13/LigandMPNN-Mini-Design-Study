import math
import csv

pdb_file = "data/target/8VEZ.pdb"
ligand_name = "GG2"
cutoff = 5.0

ligand_atoms = []
protein_atoms = []

with open(pdb_file, "r") as file:
    for line in file:

        if line.startswith("HETATM"):
            residue_name = line[17:20].strip()

            if residue_name == ligand_name:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                ligand_atoms.append((x, y, z))

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

print("Ligand atoms:", len(ligand_atoms))

nearby_residues = {}

for atom in protein_atoms:

    chain, residue_name, residue_number, atom_name, x, y, z = atom

    minimum_distance = None

    for ligand_atom in ligand_atoms:

        lx, ly, lz = ligand_atom

        distance = math.sqrt(
            (x - lx) ** 2
            + (y - ly) ** 2
            + (z - lz) ** 2
        )

        if minimum_distance is None or distance < minimum_distance:
            minimum_distance = distance

    if minimum_distance <= cutoff:

        residue_key = (
            chain,
            residue_number,
            residue_name
        )

        if (
            residue_key not in nearby_residues
            or minimum_distance < nearby_residues[residue_key]
        ):
            nearby_residues[residue_key] = minimum_distance


print(f"\nResidues within {cutoff} Å of {ligand_name}:")

rows = []

for residue, distance in sorted(
    nearby_residues.items(),
    key=lambda x: float(x[0][1])
):

    chain, residue_number, residue_name = residue

    print(
        chain,
        residue_number,
        residue_name,
        round(distance, 2),
        "Å"
    )

    rows.append({
        "chain": chain,
        "residue_number": residue_number,
        "residue_name": residue_name,
        "minimum_distance_A": round(distance, 2)
    })


output_file = "results/8VEZ_pocket_residues.csv"

with open(output_file, "w", newline="") as csvfile:

    writer = csv.DictWriter(
        csvfile,
        fieldnames=[
            "chain",
            "residue_number",
            "residue_name",
            "minimum_distance_A"
        ]
    )

    writer.writeheader()
    writer.writerows(rows)

print("\nPocket residue count:", len(rows))
print("Saved:", output_file)