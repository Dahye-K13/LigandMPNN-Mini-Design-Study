# Day 11 - Mutation pattern visualization

## Main result
Across these ten designs, the mean pocket mutation rate was
26.92%, compared with 46.52% outside the pocket.
The pocket rate was lower in 10/10 designs.
The mean paired difference was
-19.59 percentage points.

## Paired comparison
Each line connects the pocket and non-pocket rates for the same design.
Black diamonds show the means.

![Paired rates](../figures/day11_paired_mutation_rates.png)

## Residue-level patterns
The x-axis uses FASTA sequence positions, not PDB residue numbers.
Pocket residues conserved in every design: 14/26.

![Residue frequencies](../figures/day11_residue_mutation_frequency.png)

## Five highest-frequency pocket rows
All pocket residues are saved in results/day11_ranked_pocket_residues.csv.
Ties are ordered by sequence position; the top-five boundary may split ties.

| PDB residue | Sequence position | Mutated designs | Frequency |
|---|---:|---:|---:|
| A:I11 | 12 | 10/10 | 100% |
| A:F89 | 90 | 10/10 | 100% |
| A:L91 | 92 | 10/10 | 100% |
| A:S55 | 56 | 9/10 | 90% |
| A:A56 | 57 | 8/10 | 80% |

## Interpretation and limits
Pocket residues were more conserved in this set of generated sequences.
A frequently mutated position is a candidate for further inspection,
not evidence of improved binding or an experimentally validated hotspot.
These samples share one backbone and sampling setup.
No statistical significance or causal ligand-context effect is claimed.
The pocket uses the Day 10 heavy-atom distance definition (5 angstroms);
it is not assumed identical to LigandMPNN's internal scoring mask.
