pdb_file = "data/raw/1BC8.pdb"

atom_count = 0
hetatm_count = 0

with open(pdb_file, "r") as file:
    for line in file:

        if line.startswith("ATOM"):
            atom_count += 1

        elif line.startswith("HETATM"):
            hetatm_count += 1

print("ATOM count:", atom_count)
print("HETATM count:", hetatm_count)