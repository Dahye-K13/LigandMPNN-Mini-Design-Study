from pathlib import Path
import csv

# LigandMPNN output FASTA
fasta_file = Path.home() / "LigandMPNN" / "outputs" / "ligandmpnn_day3" / "seqs" / "1BC8.fa"

records = []

header = None
sequence_lines = []

with open(fasta_file, "r") as file:
    for line in file:
        line = line.strip()

        if line.startswith(">"):
            if header is not None:
                sequence = "".join(sequence_lines)
                records.append((header, sequence))

            header = line
            sequence_lines = []

        else:
            sequence_lines.append(line)

# 마지막 sequence 저장
if header is not None:
    sequence = "".join(sequence_lines)
    records.append((header, sequence))

# 첫 번째 = native, 두 번째 = designed
native_header, native = records[0]
design_header, design = records[1]

print("Native sequence length:", len(native))
print("Design sequence length:", len(design))

mutations = []
same_count = 0

for position, (wt, mut) in enumerate(zip(native, design), start=1):

    if wt == mut:
        same_count += 1

    else:
        mutations.append({
            "position": position,
            "native": wt,
            "design": mut
        })

sequence_identity = same_count / len(native)

print("Mutation count:", len(mutations))
print("Sequence identity:", round(sequence_identity, 4))
print("Sequence identity (%):", round(sequence_identity * 100, 2))

print("\nFirst 10 mutations:")

for mutation in mutations[:10]:
    print(
        mutation["position"],
        mutation["native"],
        "->",
        mutation["design"]
    )

# CSV 저장
output_file = "results/day7_mutations.csv"

with open(output_file, "w", newline="") as csvfile:
    fieldnames = ["position", "native", "design"]

    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()
    writer.writerows(mutations)

print("\nSaved:", output_file)