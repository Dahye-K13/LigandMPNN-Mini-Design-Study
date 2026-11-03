# Multiple LigandMPNN Designs

## Target
- PDB: 8VEZ
- Protein: apx1049
- Ligand: Apixaban (GG2)

## Design settings
- Temperature: 0.1
- Seed: 111
- Ligand context: enabled
- Ligand cutoff distance: 5.0 Å
- Batch size: 1
- Number of batches: 10

## Results
- Generated designs: 10
- Reference sequence length: 115 residues
- Summary: results/design_summary.csv

The summary includes overall confidence, ligand confidence,
sequence recovery, mutation count, whole-sequence identity,
and the designed sequence.

Mutation count and sequence identity were calculated against
the first (reference) sequence in the LigandMPNN FASTA.
Sequence identity is reported as a fraction from 0 to 1.

## Troubleshooting
The summary script initially exited without output because
the saved Python file was empty (0 bytes).
After saving the code and running it, the CSV was generated
with one header row and ten design rows.

## Limitation
LigandMPNN confidence scores do not directly measure
experimental protein stability or binding affinity.