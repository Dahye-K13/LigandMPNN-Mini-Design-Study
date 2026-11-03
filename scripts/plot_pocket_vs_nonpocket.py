import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


# -----------------------------
# 1. Load data
# -----------------------------

csv_file = "results/pocket_mutation_summary.csv"

df = pd.read_csv(csv_file)

print("Loaded:", csv_file)
print("Columns:", df.columns.tolist())


# -----------------------------
# 2. Find required columns
# -----------------------------

pocket_col = None
nonpocket_col = None
design_col = None


for col in df.columns:

    normalized = (
        col.lower()
        .replace("-", "_")
        .replace(" ", "_")
    )

    if (
        "pocket" in normalized
        and "mutation" in normalized
        and "rate" in normalized
    ):
        if "non" in normalized:
            nonpocket_col = col
        else:
            pocket_col = col

    if normalized in ["design_id", "design", "id"]:
        design_col = col


if pocket_col is None or nonpocket_col is None:
    raise ValueError(
        "Could not identify pocket/non-pocket mutation-rate columns."
    )


print("Pocket column:", pocket_col)
print("Non-pocket column:", nonpocket_col)


# -----------------------------
# 3. Prepare values
# -----------------------------

pocket = df[pocket_col].astype(float)
nonpocket = df[nonpocket_col].astype(float)


# Convert 0–1 values to %
if max(pocket.max(), nonpocket.max()) <= 1:
    pocket = pocket * 100
    nonpocket = nonpocket * 100


if design_col is not None:
    labels = [
        f"Design {x}"
        for x in df[design_col]
    ]
else:
    labels = [
        f"Design {i}"
        for i in range(1, len(df) + 1)
    ]


# -----------------------------
# 4. Mean mutation rates
# -----------------------------

mean_pocket = pocket.mean()
mean_nonpocket = nonpocket.mean()

print()
print("Mean pocket mutation rate:",
      round(mean_pocket, 2), "%")

print("Mean non-pocket mutation rate:",
      round(mean_nonpocket, 2), "%")


# -----------------------------
# 5. Plot
# -----------------------------

x = list(range(len(df)))
width = 0.38

plt.figure(figsize=(11, 6))

plt.bar(
    [i - width / 2 for i in x],
    pocket,
    width=width,
    label="Pocket"
)

plt.bar(
    [i + width / 2 for i in x],
    nonpocket,
    width=width,
    label="Non-pocket"
)

plt.xticks(
    x,
    labels,
    rotation=45,
    ha="right"
)

plt.xlabel("LigandMPNN design")
plt.ylabel("Mutation rate (%)")

plt.title(
    "Pocket vs Non-pocket Mutation Rates Across LigandMPNN Designs"
)

plt.legend()
plt.tight_layout()


# -----------------------------
# 6. Save
# -----------------------------

figures_dir = Path("figures")
figures_dir.mkdir(exist_ok=True)

output_file = (
    figures_dir
    / "pocket_vs_nonpocket_mutation_rates.png"
)

plt.savefig(
    output_file,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print()
print("Saved:", output_file)
