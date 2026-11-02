"""Compare mutation rates around apixaban and outside its pocket."""
import csv
from pathlib import Path
from statistics import mean

ROOT = Path(__file__).resolve().parents[1]
PDB = Path.home() / "LigandMPNN/inputs/8VEZ.pdb"
FASTA = Path.home() / "LigandMPNN/outputs/8VEZ_day4/seqs/8VEZ.fa"
CUTOFF = 5.0
AA = dict(zip(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
    "ARNDCQEGHILKMFPSTWYV"
))


def save_csv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    # Read native sequence and generated designs.
    records = []
    for line in FASTA.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            records.append([line[1:], ""])
        elif records:
            records[-1][1] += line
        else:
            raise ValueError("FASTA sequence appears before its header.")

    if len(records) != 11:
        raise ValueError("Expected one reference and ten designs.")
    native = records[0][1]
    designs = []
    for header, seq in records[1:]:
        meta = dict(x.strip().split("=", 1)
                    for x in header.split(",") if "=" in x)
        if len(seq) != len(native):
            raise ValueError("Design and reference lengths differ.")
        designs.append((meta["id"], seq))
    if len({design_id for design_id, _ in designs}) != 10:
        raise ValueError("Design IDs are not unique.")

    # Read heavy atoms from the first PDB model.
    chains = {}
    ligands = {}
    for line in PDB.read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if line[:6].strip() not in {"ATOM", "HETATM"}:
            continue
        if line[16] not in {" ", "A"}:
            continue
        atom = line[12:16].strip()
        element = line[76:78].strip()
        if element in {"H", "D"} or (
            not element and atom.lstrip("0123456789").startswith(("H", "D"))
        ):
            continue
        name = line[17:20].strip()
        chain = line[21]
        key = (int(line[22:26]), line[26].strip())
        xyz = tuple(float(line[start:start+8]) for start in (30, 38, 46))
        if name == "GG2":
            ligands.setdefault((chain, *key), []).append(xyz)
        elif line.startswith("ATOM") and name in AA:
            residues = chains.setdefault(chain, {})
            residue = residues.setdefault(
                key, {"aa": AA[name], "atoms": [], "has_ca": False}
            )
            residue["atoms"].append(xyz)
            residue["has_ca"] |= atom == "CA"

    # Require an exact, unambiguous chain-to-FASTA match.
    matches = []
    for chain, residues in chains.items():
        if not all(r["has_ca"] for r in residues.values()):
            continue
        sequence = "".join(r["aa"] for r in residues.values())
        if sequence == native:
            matches.append((chain, residues))
    if len(matches) != 1:
        raise ValueError(
            f"Found {len(matches)} exact PDB-chain matches to FASTA; "
            "manual chain/sequence mapping is required."
        )
    if len(ligands) != 1:
        raise ValueError(
            f"Found {len(ligands)} GG2 ligands; select the intended ligand first."
        )

    chain, residues = matches[0]
    ligand_key, ligand_atoms = next(iter(ligands.items()))
    positions = []
    for index, ((number, insertion), residue) in enumerate(residues.items()):
        distance = min(
            sum((a - b) ** 2 for a, b in zip(protein_atom, ligand_atom))
            for protein_atom in residue["atoms"]
            for ligand_atom in ligand_atoms
        ) ** 0.5
        positions.append({
            "sequence_position": index + 1,
            "chain": chain,
            "residue_number": number,
            "insertion_code": insertion,
            "native_aa": residue["aa"],
            "minimum_distance_A": distance,
            "region": "pocket" if distance <= CUTOFF else "non_pocket",
        })

    pocket = [i for i, row in enumerate(positions) if row["region"] == "pocket"]
    outside = [i for i, row in enumerate(positions) if row["region"] == "non_pocket"]
    if not pocket or not outside:
        raise ValueError("Both pocket and non-pocket regions must be nonempty.")

    summaries = []
    for design_id, seq in designs:
        p_mut = sum(seq[i] != native[i] for i in pocket)
        n_mut = sum(seq[i] != native[i] for i in outside)
        assert p_mut + n_mut == sum(a != b for a, b in zip(seq, native))
        summaries.append({
            "design_id": design_id,
            "pocket_residues": len(pocket),
            "pocket_mutations": p_mut,
            "pocket_mutation_rate": p_mut / len(pocket),
            "non_pocket_residues": len(outside),
            "non_pocket_mutations": n_mut,
            "non_pocket_mutation_rate": n_mut / len(outside),
            "pocket_minus_non_pocket": p_mut / len(pocket) - n_mut / len(outside),
        })

    for i, row in enumerate(positions):
        count = sum(seq[i] != native[i] for _, seq in designs)
        row["mutated_designs"] = count
        row["mutation_frequency"] = count / len(designs)

    p_mean = mean(r["pocket_mutation_rate"] for r in summaries)
    n_mean = mean(r["non_pocket_mutation_rate"] for r in summaries)
    lower_count = sum(r["pocket_minus_non_pocket"] < 0 for r in summaries)
    report = f"""# Day 10 - Pocket versus non-pocket mutations

## Inputs and mapping
- Structure: ~/LigandMPNN/inputs/8VEZ.pdb
- Designs: ~/LigandMPNN/outputs/8VEZ_day4/seqs/8VEZ.fa
- Designs analyzed: {len(designs)}
- Matched protein chain: {chain}
- Reference length: {len(native)}
- Ligand identifier (chain, residue number, insertion code): {ligand_key}
- PDB chain sequence exactly matches the FASTA reference.
- First PDB model; blank/A alternate locations; standard ATOM residues.

## Pocket definition
A residue is in the pocket when any of its heavy atoms is within
{CUTOFF} angstroms of any GG2 heavy atom.
This definition was recalculated from the PDB, independently of the old
pocket CSV and LigandMPNN's internal ligand-confidence residue selection.

## Results
- Pocket residues: {len(pocket)}
- Non-pocket residues: {len(outside)}
- Mean pocket mutation rate: {p_mean:.2%}
- Mean non-pocket mutation rate: {n_mean:.2%}
- Mean pocket minus non-pocket rate: {(p_mean - n_mean) * 100:.2f} percentage points
- Designs with a lower pocket mutation rate: {lower_count}/{len(designs)}

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
"""
    save_csv(ROOT / "results/day10_pocket_mutation_summary.csv", summaries)
    save_csv(ROOT / "results/day10_residue_mutation_frequency.csv", positions)
    (ROOT / "docs/day10_pocket_mutations.md").write_text(report)
    print(report)


if __name__ == "__main__":
    main()
