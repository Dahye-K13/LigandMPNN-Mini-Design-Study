"""Visualize pocket and non-pocket mutation patterns."""

import csv
from pathlib import Path
from statistics import mean

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter


ROOT = Path(__file__).resolve().parents[1]


def read_csv(name):
    """Read a CSV file from the results directory."""

    with (ROOT / "results" / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def main():

    # --------------------------------------------------
    # 1. Load mutation analysis results
    # --------------------------------------------------

    summary = read_csv(
        "pocket_mutation_summary.csv"
    )

    residues = read_csv(
        "residue_mutation_frequency.csv"
    )


    if len(summary) != 10:
        raise ValueError(
            "Expected 10 generated designs."
        )


    if len(residues) != 115:
        raise ValueError(
            "Expected 115 residue rows."
        )


    # --------------------------------------------------
    # 2. Extract mutation rates
    # --------------------------------------------------

    pocket = [
        float(row["pocket_mutation_rate"])
        for row in summary
    ]

    outside = [
        float(row["non_pocket_mutation_rate"])
        for row in summary
    ]


    positions = [
        int(row["sequence_position"])
        for row in residues
    ]

    frequencies = [
        float(row["mutation_frequency"])
        for row in residues
    ]


    if positions != list(range(1, 116)):
        raise ValueError(
            "Residue rows must be in sequence order."
        )


    if not all(
        0 <= value <= 1
        for value
        in pocket + outside + frequencies
    ):
        raise ValueError(
            "Mutation rates must be between 0 and 1."
        )


    # --------------------------------------------------
    # 3. Prepare figure directory
    # --------------------------------------------------

    figures = ROOT / "figures"

    figures.mkdir(
        exist_ok=True
    )


    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


    # --------------------------------------------------
    # 4. Paired pocket vs non-pocket mutation rates
    # --------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(6, 5)
    )


    for pocket_rate, non_pocket_rate in zip(
        pocket,
        outside
    ):

        ax.plot(
            [0, 1],
            [
                pocket_rate,
                non_pocket_rate
            ],
            color="#94a3b8",
            alpha=0.65,
            lw=1
        )

        ax.scatter(
            [0, 1],
            [
                pocket_rate,
                non_pocket_rate
            ],
            c=[
                "#0f766e",
                "#c2410c"
            ],
            s=30
        )


    ax.scatter(
        [0, 1],
        [
            mean(pocket),
            mean(outside)
        ],
        marker="D",
        s=110,
        c="black",
        label="Mean",
        zorder=5
    )


    ax.set_xticks(
        [0, 1],
        [
            "Pocket\n26 residues",
            "Non-pocket\n89 residues"
        ]
    )

    ax.set_xlim(
        -0.35,
        1.35
    )

    ax.set_ylim(
        0,
        1
    )

    ax.yaxis.set_major_formatter(
        PercentFormatter(1)
    )

    ax.set_ylabel(
        "Mutation rate within each region"
    )

    ax.set_title(
        "8VEZ: paired mutation rates across 10 designs"
    )

    ax.legend(
        frameon=False
    )

    ax.grid(
        axis="y",
        alpha=0.2
    )

    fig.tight_layout()


    paired_figure = (
        figures
        / "paired_mutation_rates.png"
    )

    fig.savefig(
        paired_figure,
        dpi=200
    )

    plt.close(fig)


    # --------------------------------------------------
    # 5. Mutation frequency by residue
    # --------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(12, 4)
    )


    colors = [
        "#0f766e"
        if row["region"] == "pocket"
        else "#cbd5e1"
        for row in residues
    ]


    ax.bar(
        positions,
        frequencies,
        color=colors,
        width=0.9
    )


    # Dummy bars for legend
    ax.bar(
        [],
        [],
        color="#0f766e",
        label="Pocket"
    )

    ax.bar(
        [],
        [],
        color="#cbd5e1",
        label="Non-pocket"
    )


    ax.set_xlim(
        0,
        116
    )

    ax.set_ylim(
        0,
        1.05
    )

    ax.yaxis.set_major_formatter(
        PercentFormatter(1)
    )

    ax.set_xlabel(
        "Sequence position (1-based; not PDB residue number)"
    )

    ax.set_ylabel(
        "Fraction of designs mutated"
    )

    ax.set_title(
        "8VEZ: mutation frequency at each residue"
    )

    ax.legend(
        frameon=False
    )

    fig.tight_layout()


    frequency_figure = (
        figures
        / "residue_mutation_frequency.png"
    )

    fig.savefig(
        frequency_figure,
        dpi=200
    )

    plt.close(fig)


    # --------------------------------------------------
    # 6. Rank pocket residues by mutation frequency
    # --------------------------------------------------

    ranked = sorted(
        (
            row
            for row in residues
            if row["region"] == "pocket"
        ),
        key=lambda row: (
            -float(
                row["mutation_frequency"]
            ),
            int(
                row["sequence_position"]
            )
        )
    )


    if len(ranked) != 26:
        raise ValueError(
            "Expected 26 pocket residues."
        )


    ranked_output = (
        ROOT
        / "results/ranked_pocket_residues.csv"
    )


    with ranked_output.open(
        "w",
        newline=""
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=list(
                ranked[0]
            )
        )

        writer.writeheader()
        writer.writerows(
            ranked
        )


    # --------------------------------------------------
    # 7. Prepare top-five residue table
    # --------------------------------------------------

    top_rows = []


    for row in ranked[:5]:

        label = (
            f'{row["chain"]}:'
            f'{row["native_aa"]}'
            f'{row["residue_number"]}'
            f'{row["insertion_code"]}'
        )


        top_rows.append(
            f'| {label} | '
            f'{row["sequence_position"]} | '
            f'{int(row["mutated_designs"])}/10 | '
            f'{float(row["mutation_frequency"]):.0%} |'
        )


    # --------------------------------------------------
    # 8. Summary statistics
    # --------------------------------------------------

    lower_count = sum(
        p < n
        for p, n
        in zip(
            pocket,
            outside
        )
    )


    conserved_count = sum(
        int(
            row["mutated_designs"]
        ) == 0
        for row in ranked
    )


    paired_difference = mean(
        p - n
        for p, n
        in zip(
            pocket,
            outside
        )
    )


    # --------------------------------------------------
    # 9. Create report
    # --------------------------------------------------

    report = f"""# Mutation pattern visualization

## Main result

Across these ten designs, the mean pocket mutation rate was
{mean(pocket):.2%}, compared with {mean(outside):.2%} outside the pocket.

The pocket mutation rate was lower in {lower_count}/10 designs.

The mean paired difference was
{100 * paired_difference:.2f} percentage points.

## Paired comparison

Each line connects the pocket and non-pocket mutation rates
for the same LigandMPNN design.

Black diamonds indicate the mean values.

![Paired mutation rates](../figures/paired_mutation_rates.png)

## Residue-level patterns

The x-axis uses FASTA sequence positions rather than PDB residue numbers.

Pocket residues conserved in every design: {conserved_count}/26.

![Residue mutation frequencies](../figures/residue_mutation_frequency.png)

## Five highest-frequency pocket residues

All pocket residues ranked by mutation frequency are saved in:

results/ranked_pocket_residues.csv

Ties are ordered by sequence position.

| PDB residue | Sequence position | Mutated designs | Frequency |
|---|---:|---:|---:|
{chr(10).join(top_rows)}

## Interpretation and limitations

Pocket residues were more conserved than non-pocket residues
across this set of generated sequences.

However, individual pocket positions showed substantially different
mutation frequencies.

A frequently mutated residue should be interpreted as a candidate
for further structural inspection rather than evidence of improved
binding or an experimentally validated binding hotspot.

All designs were generated from a single protein backbone using
one sampling configuration.

No statistical significance, experimental binding affinity,
protein stability, or causal effect of ligand context is claimed.

The ligand-binding pocket was defined using a 5 Å heavy-atom
distance cutoff calculated directly from the 8VEZ structure.
"""


    report_output = (
        ROOT
        / "docs/mutation_visualization.md"
    )

    report_output.write_text(
        report
    )


    # --------------------------------------------------
    # 10. Terminal summary
    # --------------------------------------------------

    print(
        f"Mean pocket rate: "
        f"{mean(pocket):.2%}"
    )

    print(
        f"Mean non-pocket rate: "
        f"{mean(outside):.2%}"
    )

    print(
        "Pocket residues unchanged "
        f"in all designs: "
        f"{conserved_count}/26"
    )

    print(
        "\nHighest-frequency pocket residues:"
    )

    print(
        "\n".join(
            top_rows
        )
    )

    print(
        "\nSaved:"
    )

    print(
        f"- {paired_figure}"
    )

    print(
        f"- {frequency_figure}"
    )

    print(
        f"- {ranked_output}"
    )

    print(
        f"- {report_output}"
    )


if __name__ == "__main__":
    main()