#!/usr/bin/env python3

"""
Consolidate preliminary Sanger SNV and indel observations.

Preserves one row per raw read-level observation.
Does not normalize indels, infer genotypes, or validate variants.
"""

import argparse
import csv
from collections import Counter
from pathlib import Path


COLUMNS = [
    "Sample", "Gene", "Target", "Read_ID", "Direction",
    "Chromosome", "Position_GRCh38_1based",
    "REF", "ALT", "Variant_type",
    "Observation_quality",
    "Evidence_flag", "Validation_status",
    "Evidence_supporting_read_IDs",
    "Evidence_opposing_read_IDs",
    "Notes",
]


def load_csv(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def allele_key(row):
    return (
        row["Sample"],
        row["Gene"],
        row["Target"],
        row["Chromosome"],
        str(int(row["Position_GRCh38_1based"])),
        row["REF"].upper(),
        row["ALT"].upper(),
    )


def indel_key(row):
    return (row["Read_ID"],) + allele_key(row)


def write_csv(path, columns, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Consolidate unvalidated Sanger observations."
    )
    parser.add_argument("--raw-snvs", required=True, type=Path)
    parser.add_argument("--snv-evidence", required=True, type=Path)
    parser.add_argument("--raw-indels", required=True, type=Path)
    parser.add_argument("--indel-evidence", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    paths = {
        name: path.expanduser().resolve()
        for name, path in {
            "raw_snvs": args.raw_snvs,
            "snv_evidence": args.snv_evidence,
            "raw_indels": args.raw_indels,
            "indel_evidence": args.indel_evidence,
        }.items()
    }

    for name, path in paths.items():
        if not path.is_file():
            parser.error(f"Missing {name}: {path}")

    datasets = {
        name: [
            row for row in load_csv(path)
            if row["Sample"] == args.sample
        ]
        for name, path in paths.items()
    }

    raw_snvs = datasets["raw_snvs"]
    snv_evidence = datasets["snv_evidence"]
    raw_indels = datasets["raw_indels"]
    indel_evidence = datasets["indel_evidence"]

    snv_lookup = {}
    for row in snv_evidence:
        key = allele_key(row)
        if key in snv_lookup:
            raise ValueError(f"Duplicate SNV evidence: {key}")
        snv_lookup[key] = row

    indel_lookup = {}
    for row in indel_evidence:
        key = indel_key(row)
        if key in indel_lookup:
            raise ValueError(f"Duplicate indel evidence: {key}")
        indel_lookup[key] = row

    candidates = []
    used_snv_keys = set()
    used_indel_keys = set()

    for row in raw_snvs:
        key = allele_key(row)
        if key not in snv_lookup:
            raise ValueError(f"Missing SNV evidence: {key}")

        evidence = snv_lookup[key]
        used_snv_keys.add(key)

        candidates.append({
            "Sample": row["Sample"],
            "Gene": row["Gene"],
            "Target": row["Target"],
            "Read_ID": row["Read_ID"],
            "Direction": row["Direction"],
            "Chromosome": row["Chromosome"],
            "Position_GRCh38_1based": row["Position_GRCh38_1based"],
            "REF": row["REF"],
            "ALT": row["ALT"],
            "Variant_type": "SNV",
            "Observation_quality": row["Base_quality"],
            "Evidence_flag": evidence["Evidence_flag"],
            "Validation_status": "UNVALIDATED",
            "Evidence_supporting_read_IDs": evidence["Supporting_read_IDs"],
            "Evidence_opposing_read_IDs": evidence["Opposing_read_IDs"],
            "Notes": (
                "Q20 evidence evaluated across passing reads; "
                "no genotype inferred"
            ),
        })

    for row in raw_indels:
        key = indel_key(row)
        if key not in indel_lookup:
            raise ValueError(f"Missing indel evidence: {key}")

        evidence = indel_lookup[key]
        used_indel_keys.add(key)

        candidates.append({
            "Sample": row["Sample"],
            "Gene": row["Gene"],
            "Target": row["Target"],
            "Read_ID": row["Read_ID"],
            "Direction": row["Direction"],
            "Chromosome": row["Chromosome"],
            "Position_GRCh38_1based": row["Position_GRCh38_1based"],
            "REF": row["REF"],
            "ALT": row["ALT"],
            "Variant_type": row["Variant_type"],
            "Observation_quality": (
                f"left={row['Left_flank_quality']};"
                f"right={row['Right_flank_quality']}"
            ),
            "Evidence_flag": evidence["Evidence_flag"],
            "Validation_status": "UNVALIDATED",
            "Evidence_supporting_read_IDs": evidence[
                "Exact_deletion_support_reads"
            ],
            "Evidence_opposing_read_IDs": evidence[
                "Intact_reference_Q20_reads"
            ],
            "Notes": (
                "Indel not normalized; flanking base quality "
                "does not establish independent confirmation"
            ),
        })

    if used_snv_keys != set(snv_lookup):
        raise ValueError("SNV evidence contains unmatched alleles")

    if used_indel_keys != set(indel_lookup):
        raise ValueError("Indel evidence contains unmatched observations")

    if len(candidates) != len(raw_snvs) + len(raw_indels):
        raise RuntimeError("Candidate count mismatch")

    candidates.sort(
        key=lambda row: (
            row["Chromosome"],
            int(row["Position_GRCh38_1based"]),
            row["Variant_type"],
            row["Read_ID"],
        )
    )

    summary_counts = Counter(
        (row["Variant_type"], row["Evidence_flag"])
        for row in candidates
    )

    summary = [
        {
            "Sample": args.sample,
            "Variant_type": variant_type,
            "Evidence_flag": evidence_flag,
            "Observation_count": count,
            "Validation_status": "UNVALIDATED",
        }
        for (variant_type, evidence_flag), count
        in sorted(summary_counts.items())
    ]

    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    candidate_path = output_dir / (
        f"{args.sample}_Preliminary_Candidates.csv"
    )
    summary_path = output_dir / (
        f"{args.sample}_Candidate_Evidence_Summary.csv"
    )

    write_csv(candidate_path, COLUMNS, candidates)
    write_csv(
        summary_path,
        [
            "Sample", "Variant_type", "Evidence_flag",
            "Observation_count", "Validation_status"
        ],
        summary,
    )

    print("\nPRELIMINARY SANGER CANDIDATE CONSOLIDATION")
    print("=" * 65)
    print("Sample:", args.sample)
    print("Raw SNV observations:", len(raw_snvs))
    print("Distinct SNV alleles:", len(used_snv_keys))
    print("Raw indel observations:", len(raw_indels))
    print("Total candidate observation rows:", len(candidates))

    print("\nEvidence summary:")
    for row in summary:
        print(
            f"  {row['Variant_type']} | "
            f"{row['Evidence_flag']}: "
            f"{row['Observation_count']}"
        )

    print("\nSaved:")
    print(" ", candidate_path)
    print(" ", summary_path)

    print("\nAll candidates remain UNVALIDATED.")
    print("Rows represent observations, not distinct confirmed variants.")


if __name__ == "__main__":
    main()
