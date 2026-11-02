"""Summarize LigandMPNN designs relative to the first FASTA record."""
import csv
from pathlib import Path

def main():
    source = Path.home() / "LigandMPNN/outputs/8VEZ_day4/seqs/8VEZ.fa"
    output = Path(__file__).resolve().parents[1] / "results/day4_design_summary.csv"
    records = []
    for line in source.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            records.append([line[1:], ""])
        elif records:
            records[-1][1] += line
        else:
            raise ValueError("Sequence before FASTA header")
    if len(records) < 2 or not records[0][1]:
        raise ValueError("Expected a reference and at least one design")
    reference = records[0][1]
    rows = []
    seen = set()
    for header, sequence in records[1:]:
        metadata = dict(part.strip().split("=", 1) for part in header.split(",") if "=" in part)
        design_id = metadata["id"]
        if design_id in seen:
            raise ValueError(f"Duplicate design ID: {design_id}")
        seen.add(design_id)
        if len(sequence) != len(reference):
            raise ValueError(f"Design {design_id}: length differs from reference")
        mutations = sum(a != b for a, b in zip(reference, sequence))
        rows.append({
            "design_id": design_id,
            "overall_confidence": float(metadata["overall_confidence"]),
            "ligand_confidence": float(metadata["ligand_confidence"]),
            "seq_rec": float(metadata["seq_rec"]),
            "mutation_count": mutations,
            "sequence_identity": (len(reference) - mutations) / len(reference),
            "sequence": sequence,
        })
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} designs to {output}")
    print(f"Reference length: {len(reference)}; sequence_identity is a fraction (0-1).")

if __name__ == "__main__":
    main()