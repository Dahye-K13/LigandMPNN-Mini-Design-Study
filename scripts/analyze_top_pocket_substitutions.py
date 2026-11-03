import pandas as pd
from collections import Counter


# --------------------------------
# 1. Load data
# --------------------------------

ranked_file = "results/day4_ranked_pocket_residues.csv"
design_file = "results/day4_design_summary.csv"

ranked = pd.read_csv(ranked_file)
designs = pd.read_csv(design_file)


# Top 5 pocket residues
top5 = ranked.head(5).copy()


# --------------------------------
# 2. Analyze substitutions
# --------------------------------

results = []

for _, row in top5.iterrows():

    sequence_position = int(row["sequence_position"])
    residue_number = int(row["residue_number"])
    native_aa = row["native_aa"]

    designed_aas = []

    for sequence in designs["sequence"]:

        # sequence_position is 1-based
        aa = sequence[sequence_position - 1]

        designed_aas.append(aa)


    all_counts = Counter(designed_aas)

    mutant_aas = [
        aa
        for aa in designed_aas
        if aa != native_aa
    ]

    mutant_counts = Counter(mutant_aas)


    if mutant_counts:

        most_common_mutant, count = (
            mutant_counts.most_common(1)[0]
        )

    else:

        most_common_mutant = "-"
        count = 0


    results.append(
        {
            "residue": f"{native_aa}{residue_number}",
            "sequence_position": sequence_position,
            "native_aa": native_aa,
            "mutation_frequency": row["mutation_frequency"],
            "designed_aa_counts": "; ".join(
                f"{aa}:{count}"
                for aa, count in all_counts.most_common()
            ),
            "mutant_aa_counts": "; ".join(
                f"{aa}:{count}"
                for aa, count in mutant_counts.most_common()
            ),
            "most_common_mutant": most_common_mutant,
            "most_common_mutant_count": count
        }
    )


# --------------------------------
# 3. Save
# --------------------------------

result_df = pd.DataFrame(results)

output_file = (
    "results/day4_top5_pocket_substitutions.csv"
)

result_df.to_csv(
    output_file,
    index=False
)


print(
    result_df.to_string(index=False)
)

print("\nSaved:", output_file)
