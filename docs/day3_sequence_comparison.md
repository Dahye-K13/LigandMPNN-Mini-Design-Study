# Day 7 - LigandMPNN Sequence Comparison

## Input
PDB: 1BC8

## Generated Design
One sequence generated using LigandMPNN.

## Metrics
- overall_confidence:0.4784
- ligand_confidence:0.5487
- seq_rec:0.4839

## Sequence Comparison
- Total sequence length:93 residues
- Number of mutations:48
- Whole-sequence identity:48.39%

## Example Mutations
## Example Mutations

- D2K: Aspartate (negative) → Lysine (positive), representing a charge reversal.
- T6S: Threonine → Serine; both are polar residues with hydroxyl groups, representing a relatively conservative substitution.
- W8H: Tryptophan → Histidine; a bulky hydrophobic aromatic residue was replaced by a smaller polar/basic aromatic residue.
- K17D: Lysine (positive) → Aspartate (negative), representing another charge reversal.

## Interpretation
LigandMPNN redesigned 48 of 93 residues, resulting in a whole-sequence identity of 48.39% relative to the native sequence. The substitutions included both relatively conservative changes and substantial physicochemical changes such as charge reversals.

## Limitation
LigandMPNN confidence scores represent model confidence and do not directly demonstrate experimental binding affinity or protein stability. Further structural and experimental validation is required.