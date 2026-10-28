fasta_file = "data/raw/example.fasta"

with open(fasta_file, "r") as file:
    for line in file:

        line = line.strip()

        if line.startswith(">"):
            print("Header:", line)

        else:
            print("Sequence:", line)
            print("Length:", len(line))