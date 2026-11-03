"""Plot Day 10 mutation rates and rank pocket mutation frequencies."""
import csv
from pathlib import Path
from statistics import mean

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

ROOT = Path(__file__).resolve().parents[1]


def read_csv(name):
    with (ROOT / "results" / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def main():
    summary = read_csv("day10_pocket_mutation_summary.csv")
    residues = read_csv("day10_residue_mutation_frequency.csv")
    if len(summary) != 10 or len(residues) != 115:
        raise ValueError("Expected 10 designs and 115 residue rows.")

    pocket = [float(r["pocket_mutation_rate"]) for r in summary]
    outside = [float(r["non_pocket_mutation_rate"]) for r in summary]
    positions = [int(r["sequence_position"]) for r in residues]
    frequencies = [float(r["mutation_frequency"]) for r in residues]
    if positions != list(range(1, 116)):
        raise ValueError("Residue rows must be in sequence order.")
    if not all(0 <= x <= 1 for x in pocket + outside + frequencies):
        raise ValueError("Mutation rates must be between 0 and 1.")

    figures = ROOT / "figures"
    figures.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                         "axes.spines.right": False})

    # Each line connects the two regions of the SAME design.
    fig, ax = plt.subplots(figsize=(6, 5))
    for p, n in zip(pocket, outside):
        ax.plot([0, 1], [p, n], color="#94a3b8", alpha=0.65, lw=1)
        ax.scatter([0, 1], [p, n], c=["#0f766e", "#c2410c"], s=30)
    ax.scatter([0, 1], [mean(pocket), mean(outside)],
               marker="D", s=110, c="black", label="Mean", zorder=5)
    ax.set_xticks([0, 1], ["Pocket\n26 residues", "Non-pocket\n89 residues"])
    ax.set_xlim(-0.35, 1.35)
    ax.set_ylim(0, 1)
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.set_ylabel("Mutation rate within each region")
    ax.set_title("8VEZ: paired mutation rates across 10 designs")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(figures / "day11_paired_mutation_rates.png", dpi=200)
    plt.close(fig)

    # Frequency across designs at each sequence position.
    fig, ax = plt.subplots(figsize=(12, 4))
    colors = ["#0f766e" if r["region"] == "pocket" else "#cbd5e1"
              for r in residues]
    ax.bar(positions, frequencies, color=colors, width=0.9)
    ax.bar([], [], color="#0f766e", label="Pocket")
    ax.bar([], [], color="#cbd5e1", label="Non-pocket")
    ax.set_xlim(0, 116)
    ax.set_ylim(0, 1.05)
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.set_xlabel("Sequence position (1-based; not PDB residue number)")
    ax.set_ylabel("Fraction of designs mutated")
    ax.set_title("8VEZ: mutation frequency at each residue")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(figures / "day11_residue_mutation_frequency.png", dpi=200)
    plt.close(fig)

    # Save ALL pocket residues, ranked by mutation frequency.
    ranked = sorted(
        (r for r in residues if r["region"] == "pocket"),
        key=lambda r: (-float(r["mutation_frequency"]),
                       int(r["sequence_position"]))
    )
    if len(ranked) != 26:
        raise ValueError("Expected 26 pocket residues.")
    output = ROOT / "results/day11_ranked_pocket_residues.csv"
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(ranked[0]))
        writer.writeheader()
        writer.writerows(ranked)

    top_rows = []
    for r in ranked[:5]:
        label = (f'{r["chain"]}:{r["native_aa"]}'
                 f'{r["residue_number"]}{r["insertion_code"]}')
        top_rows.append(
            f'| {label} | {r["sequence_position"]} | '
            f'{int(r["mutated_designs"])}/10 | '
            f'{float(r["mutation_frequency"]):.0%} |'
        )

    lower = sum(p < n for p, n in zip(pocket, outside))
    conserved = sum(int(r["mutated_designs"]) == 0 for r in ranked)
    report = f"""# Day 11 - Mutation pattern visualization

## Main result
Across these ten designs, the mean pocket mutation rate was
{mean(pocket):.2%}, compared with {mean(outside):.2%} outside the pocket.
The pocket rate was lower in {lower}/10 designs.
The mean paired difference was
{100 * mean(p - n for p, n in zip(pocket, outside)):.2f} percentage points.

## Paired comparison
Each line connects the pocket and non-pocket rates for the same design.
Black diamonds show the means.

![Paired rates](../figures/day11_paired_mutation_rates.png)

## Residue-level patterns
The x-axis uses FASTA sequence positions, not PDB residue numbers.
Pocket residues conserved in every design: {conserved}/26.

![Residue frequencies](../figures/day11_residue_mutation_frequency.png)

## Five highest-frequency pocket rows
All pocket residues are saved in results/day11_ranked_pocket_residues.csv.
Ties are ordered by sequence position; the top-five boundary may split ties.

| PDB residue | Sequence position | Mutated designs | Frequency |
|---|---:|---:|---:|
{chr(10).join(top_rows)}

## Interpretation and limits
Pocket residues were more conserved in this set of generated sequences.
A frequently mutated position is a candidate for further inspection,
not evidence of improved binding or an experimentally validated hotspot.
These samples share one backbone and sampling setup.
No statistical significance or causal ligand-context effect is claimed.
The pocket uses the Day 10 heavy-atom distance definition (5 angstroms);
it is not assumed identical to LigandMPNN's internal scoring mask.
"""
    (ROOT / "docs/day11_mutation_visualization.md").write_text(report)
    print(f"Mean pocket rate: {mean(pocket):.2%}")
    print(f"Mean non-pocket rate: {mean(outside):.2%}")
    print(f"Pocket residues unchanged in all designs: {conserved}/26")
    print("\nHighest-frequency pocket rows:")
    print("\n".join(top_rows))
    print("\nSaved two figures, one ranked CSV, and one report.")


if __name__ == "__main__":
    main()
