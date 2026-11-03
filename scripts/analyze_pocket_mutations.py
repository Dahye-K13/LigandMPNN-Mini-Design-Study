"""Compare mutation rates around apixaban and outside its pocket."""

import csv
from pathlib import Path
from statistics import mean


ROOT = Path(__file__).resolve().parents[1]

PDB = Path.home() / "LigandMPNN/inputs/8VEZ.pdb"
FASTA = Path.home() / "LigandMPNN/outputs/8VEZ_day4/seqs/8VEZ.fa"

CUTOFF = 5.0

AA = dict(
    zip(
        "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
        "ARNDCQEGHILKMFPSTWYV"
    )
)


def save_csv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0])
        )

        writer.writeheader()
        writer.writerows(rows)


def main():

    # --------------------------------------------------
    # 1. Read native sequence and generated designs
    # --------------------------------------------------

    records = []

    for line in FASTA.read_text().splitlines():

        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):

            records.append(
                [line[1:], ""]
            )

        elif records:

            records[-1][1] += line

        else:

            raise ValueError(
                "FASTA sequence appears before its header."
            )


    if len(records) != 11:

        raise ValueError(
            "Expected one reference and ten designs."
        )


    native = records[0][1]

    designs = []


    for header, seq in records[1:]:

        meta = dict(
            x.strip().split("=", 1)
            for x in header.split(",")
            if "=" in x
        )

        if len(seq) != len(native):

            raise ValueError(
                "Design and reference lengths differ."
            )

        designs.append(
            (
                meta["id"],
                seq
            )
        )


    if len(
        {
            design_id
            for design_id, _ in designs
        }
    ) != 10:

        raise ValueError(
            "Design IDs are not unique."
        )


    # --------------------------------------------------
    # 2. Read heavy atoms from the first PDB model
    # --------------------------------------------------

    chains = {}
    ligands = {}


    for line in PDB.read_text().splitlines():

        if line.startswith("ENDMDL"):
            break

        if line[:6].strip() not in {
            "ATOM",
            "HETATM"
        }:
            continue

        if line[16] not in {
            " ",
            "A"
        }:
            continue


        atom = line[12:16].strip()
        element = line[76:78].strip()


        if (
            element in {"H", "D"}
            or (
                not element
                and atom.lstrip(
                    "0123456789"
                ).startswith(
                    ("H", "D")
                )
            )
        ):
            continue


        name = line[17:20].strip()
        chain = line[21]

        key = (
            int(line[22:26]),
            line[26].strip()
        )

        xyz = tuple(
            float(
                line[start:start + 8]
            )
            for start in (
                30,
                38,
                46
            )
        )


        if name == "GG2":

            ligands.setdefault(
                (
                    chain,
                    *key
                ),
                []
            ).append(xyz)


        elif (
            line.startswith("ATOM")
            and name in AA
        ):

            residues = chains.setdefault(
                chain,
                {}
            )

            residue = residues.setdefault(
                key,
                {
                    "aa": AA[name],
                    "atoms": [],
                    "has_ca": False
                }
            )

            residue["atoms"].append(xyz)

            residue["has_ca"] |= (
                atom == "CA"
            )


    # --------------------------------------------------
    # 3. Match PDB chain to FASTA reference
    # --------------------------------------------------

    matches = []


    for chain, residues in chains.items():

        if not all(
            residue["has_ca"]
            for residue in residues.values()
        ):
            continue


        sequence = "".join(
            residue["aa"]
            for residue in residues.values()
        )


        if sequence == native:

            matches.append(
                (
                    chain,
                    residues
                )
            )


    if len(matches) != 1:

        raise ValueError(
            f"Found {len(matches)} exact PDB-chain matches to FASTA; "
            "manual chain/sequence mapping is required."
        )


    if len(ligands) != 1:

        raise ValueError(
            f"Found {len(ligands)} GG2 ligands; "
            "select the intended ligand first."
        )


    chain, residues = matches[0]

    ligand_key, ligand_atoms = next(
        iter(
            ligands.items()
        )
    )


    # --------------------------------------------------
    # 4. Calculate ligand distances and define pocket
    # --------------------------------------------------

    positions = []


    for index, (
        (number, insertion),
        residue
    ) in enumerate(
        residues.items()
    ):

        distance = min(

            sum(
                (a - b) ** 2
                for a, b in zip(
                    protein_atom,
                    ligand_atom
                )
            )

            for protein_atom in residue["atoms"]
            for ligand_atom in ligand_atoms

        ) ** 0.5


        positions.append(
            {
                "sequence_position": index + 1,
                "chain": chain,
                "residue_number": number,
                "insertion_code": insertion,
                "native_aa": residue["aa"],
                "minimum_distance_A": distance,
                "region": (
                    "pocket"
                    if distance <= CUTOFF
                    else "non_pocket"
                ),
            }
        )


    pocket = [
        i
        for i, row in enumerate(positions)
        if row["region"] == "pocket"
    ]


    outside = [
        i
        for i, row in enumerate(positions)
        if row["region"] == "non_pocket"
    ]


    if not pocket or not outside:

        raise ValueError(
            "Both pocket and non-pocket regions must be nonempty."
        )


    # --------------------------------------------------
    # 5. Calculate mutation rates for each design
    # --------------------------------------------------

    summaries = []


    for design_id, seq in designs:

        pocket_mutations = sum(
            seq[i] != native[i]
            for i in pocket
        )

        non_pocket_mutations = sum(
            seq[i] != native[i]
            for i in outside
        )


        assert (
            pocket_mutations
            + non_pocket_mutations
            ==
            sum(
                designed_aa != native_aa
                for designed_aa, native_aa
                in zip(
                    seq,
                    native
                )
            )
        )


        pocket_rate = (
            pocket_mutations
            / len(pocket)
        )

        non_pocket_rate = (
            non_pocket_mutations
            / len(outside)
        )


        summaries.append(
            {
                "design_id": design_id,
                "pocket_residues": len(pocket),
                "pocket_mutations": pocket_mutations,
                "pocket_mutation_rate": pocket_rate,
                "non_pocket_residues": len(outside),
                "non_pocket_mutations": non_pocket_mutations,
                "non_pocket_mutation_rate": non_pocket_rate,
                "pocket_minus_non_pocket": (
                    pocket_rate
                    - non_pocket_rate
                ),
            }
        )


    # --------------------------------------------------
    # 6. Calculate residue-level mutation frequency
    # --------------------------------------------------

    for i, row in enumerate(positions):

        count = sum(
            seq[i] != native[i]
            for _, seq in designs
        )

        row["mutated_designs"] = count

        row["mutation_frequency"] = (
            count
            / len(designs)
        )


    # --------------------------------------------------
    # 7. Summary statistics
    # --------------------------------------------------

    pocket_mean = mean(
        row["pocket_mutation_rate"]
        for row in summaries
    )

    non_pocket_mean = mean(
        row["non_pocket_mutation_rate"]
        for row in summaries
    )


    lower_count = sum(
        row["pocket_minus_non_pocket"] < 0
        for row in summaries
    )


    # --------------------------------------------------
    # 8. Create analysis report
    # --------------------------------------------------

    report = f"""# Pocket versus non-pocket mutations

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

A residue is classified as a pocket residue when any of its heavy atoms
is within {CUTOFF} Å of any GG2 heavy atom.

This definition was recalculated directly from the PDB structure,
independently of LigandMPNN's internal ligand-confidence residue
selection.

## Results

- Pocket residues: {len(pocket)}
- Non-pocket residues: {len(outside)}
- Mean pocket mutation rate: {pocket_mean:.2%}
- Mean non-pocket mutation rate: {non_pocket_mean:.2%}
- Mean pocket minus non-pocket rate: {(pocket_mean - non_pocket_mean) * 100:.2f} percentage points
- Designs with a lower pocket mutation rate: {lower_count}/{len(designs)}

Mutation rates use each region's own residue count as the denominator.

Per-residue mutation frequency represents the fraction of generated
designs in which a given position differs from the native sequence.

## Outputs

- results/pocket_mutation_summary.csv
- results/residue_mutation_frequency.csv

## Scope and limitations

These are descriptive results from ten generated sequences using one
protein backbone and one sampling setup.

The analysis does not establish experimental binding affinity,
protein stability, statistical significance, or a causal effect of
ligand context.

Model confidence is not analyzed in this step.
"""


    # --------------------------------------------------
    # 9. Save outputs
    # --------------------------------------------------

    save_csv(
        ROOT
        / "results/pocket_mutation_summary.csv",
        summaries
    )

    save_csv(
        ROOT
        / "results/residue_mutation_frequency.csv",
        positions
    )

    (
        ROOT
        / "docs/pocket_mutations.md"
    ).write_text(report)


    print(report)


if __name__ == "__main__":
    main()