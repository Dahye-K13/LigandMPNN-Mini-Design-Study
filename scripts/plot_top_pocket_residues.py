import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


# -----------------------------
# 1. Load ranked pocket residues
# -----------------------------

csv_file = "results/ranked_pocket_residues.csv"

df = pd.read_csv(csv_file)

# Top 10 most frequently mutated pocket residues
top = df.head(10).copy()


# -----------------------------
# 2. Convert frequency to %
# -----------------------------

top["mutation_frequency_percent"] = (
    top["mutation_frequency"] * 100
)


# -----------------------------
# 3. Create residue labels
# -----------------------------

top["label"] = (
    top["native_aa"]
    + top["residue_number"].astype(str)
)


# Reverse order for horizontal plot
top = top.iloc[::-1]


# -----------------------------
# 4. Plot
# -----------------------------

fig, ax = plt.subplots(figsize=(8, 6))

bars = ax.barh(
    top["label"],
    top["mutation_frequency_percent"]
)

ax.set_xlabel("Mutation frequency (%)")
ax.set_ylabel("Pocket residue")

ax.set_title(
    "Most Frequently Redesigned Apixaban-Pocket Residues"
)

ax.set_xlim(0, 105)


# Add values next to bars
ax.bar_label(
    bars,
    fmt="%.0f%%",
    padding=3
)

plt.tight_layout()


# -----------------------------
# 5. Save
# -----------------------------

figures_dir = Path("figures")
figures_dir.mkdir(exist_ok=True)

output_file = (
    figures_dir
    / "top_pocket_residue_mutation_frequency.png"
)

plt.savefig(
    output_file,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Saved:", output_file)
