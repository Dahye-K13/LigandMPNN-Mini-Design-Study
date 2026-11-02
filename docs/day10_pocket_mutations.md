# Day 10 - Pocket versus non-pocket mutations

## Inputs and mapping
- Structure: ~/LigandMPNN/inputs/8VEZ.pdb
- Designs: ~/LigandMPNN/outputs/8VEZ_day4/seqs/8VEZ.fa
- Designs analyzed: 10
- Matched protein chain: A
- Reference length: 115
- Ligand identifier (chain, residue number, insertion code): ('A', 201, '')
- PDB chain sequence exactly matches the FASTA reference.
- First PDB model; blank/A alternate locations; standard ATOM residues.

## Pocket definition
A residue is in the pocket when any of its heavy atoms is within
5.0 angstroms of any GG2 heavy atom.
This definition was recalculated from the PDB, independently of the old
pocket CSV and LigandMPNN's internal ligand-confidence residue selection.

## Results
- Pocket residues: 26
- Non-pocket residues: 89
- Mean pocket mutation rate: 26.92%
- Mean non-pocket mutation rate: 46.52%
- Mean pocket minus non-pocket rate: -19.59 percentage points
- Designs with a lower pocket mutation rate: 10/10

Mutation rates use each region's own residue count as denominator.
Per-residue mutation frequency is the fraction of designs mutated there.

## Outputs
- results/day10_pocket_mutation_summary.csv
- results/day10_residue_mutation_frequency.csv

## Scope and limitations
These are descriptive results for ten samples from one backbone and
one sampling setup. They do not establish experimental binding affinity,
stability, statistical significance, or a causal effect of ligand context.
Model confidence is not analyzed in this step.
