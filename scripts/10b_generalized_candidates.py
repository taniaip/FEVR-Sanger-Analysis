#!/usr/bin/env python3
"""Consolidate unvalidated, read-level Sanger observations and evidence."""

import argparse
import csv
from collections import Counter
from pathlib import Path

COLUMNS = [
    "Sample", "Gene", "Target", "Read_ID", "Direction",
    "Chromosome", "Position_GRCh38_1based", "REF", "ALT",
    "Variant_type", "Observation_quality",
    "Evidence_flag", "Validation_status",
    "Evidence_supporting_read_IDs",
    "Evidence_opposing_read_IDs",
    "Parse_status", "Alignment_review_flags",
    "Indel_size_category", "Notes",
]


def load_csv(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def allele_key(row):
    return (
        row["Sample"], row["Gene"], row["Target"],
        row["Chromosome"], str(int(row["Position_GRCh38_1based"])),
        row["REF"].upper(), row["ALT"].upper(),
    )


def observation_key(row):
    return (row["Read_ID"], row["Variant_type"]) + allele_key(row)


def unique_lookup(rows, key_function, label):
    lookup = {}
    for row in rows:
        key = key_function(row)
        if key in lookup:
            raise ValueError(f"Duplicate {label}: {key}")
        lookup[key] = row
    return lookup


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-snvs", required=True, type=Path)
    parser.add_argument("--snv-evidence", required=True, type=Path)
    parser.add_argument("--raw-indels", required=True, type=Path)
    parser.add_argument("--deletion-evidence", required=True, type=Path)
    parser.add_argument("--insertion-evidence", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    paths = {
        name: path.expanduser().resolve()
        for name, path in {
            "raw_snvs": args.raw_snvs,
            "snv_evidence": args.snv_evidence,
            "raw_indels": args.raw_indels,
            "deletion_evidence": args.deletion_evidence,
            "insertion_evidence": args.insertion_evidence,
        }.items()
    }

    for name, path in paths.items():
        if not path.is_file():
            parser.error(f"Missing {name}: {path}")

    datasets = {}
    for name, path in paths.items():
        rows = load_csv(path)
        if any(row["Sample"] != args.sample for row in rows):
            raise ValueError(f"{name} contains another sample")
        datasets[name] = rows

    snvs = datasets["raw_snvs"]
    indels = datasets["raw_indels"]

    snv_lookup = unique_lookup(
        datasets["snv_evidence"], allele_key, "SNV evidence"
    )
    deletion_lookup = unique_lookup(
        datasets["deletion_evidence"], observation_key, "deletion evidence"
    )
    insertion_lookup = unique_lookup(
        datasets["insertion_evidence"], observation_key, "insertion evidence"
    )

    candidates = []
    used_snvs = set()
    used_deletions = set()
    used_insertions = set()

    for row in snvs:
        key = allele_key(row)
        if key not in snv_lookup:
            raise ValueError(f"Missing SNV evidence: {key}")

        evidence = snv_lookup[key]
        used_snvs.add(key)

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
            "Parse_status": row["Parse_status"],
            "Alignment_review_flags": row["Alignment_review_flags"],
            "Indel_size_category": row["Indel_size_category"],
            "Notes": "Q20 read-level SNV evidence; no genotype inferred",
        })

    for row in indels:
        key = observation_key(row)
        kind = row["Variant_type"]

        if kind == "DELETION":
            lookup = deletion_lookup
            used = used_deletions
            support_column = "Exact_deletion_support_reads"
            opposing_column = "Intact_reference_Q20_reads"
        elif kind == "INSERTION":
            lookup = insertion_lookup
            used = used_insertions
            support_column = "Exact_insertion_support_reads"
            opposing_column = "No_insertion_Q20_reads"
        else:
            raise ValueError(f"Unexpected variant type: {kind}")

        if key not in lookup:
            raise ValueError(f"Missing {kind} evidence: {key}")

        evidence = lookup[key]
        used.add(key)

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
            "Variant_type": kind,
            "Observation_quality": (
                f"left={row['Left_flank_quality']};"
                f"right={row['Right_flank_quality']}"
            ),
            "Evidence_flag": evidence["Evidence_flag"],
            "Validation_status": "UNVALIDATED",
            "Evidence_supporting_read_IDs": evidence[support_column],
            "Evidence_opposing_read_IDs": evidence[opposing_column],
            "Parse_status": row["Parse_status"],
            "Alignment_review_flags": row["Alignment_review_flags"],
            "Indel_size_category": row["Indel_size_category"],
            "Notes": (
                "Unnormalized indel; repeat-placement ambiguity possible; "
                "flank quality and read support do not confirm a variant"
            ),
        })

    for label, used, lookup in [
        ("SNV", used_snvs, snv_lookup),
        ("deletion", used_deletions, deletion_lookup),
        ("insertion", used_insertions, insertion_lookup),
    ]:
        if used != set(lookup):
            raise ValueError(f"Unmatched {label} evidence records")

    if len(candidates) != len(snvs) + len(indels):
        raise RuntimeError("Candidate observation count mismatch")

    candidates.sort(key=lambda row: (
        row["Chromosome"],
        int(row["Position_GRCh38_1based"]),
        row["Variant_type"],
        row["Read_ID"],
    ))

    counts = Counter(
        (row["Variant_type"], row["Evidence_flag"])
        for row in candidates
    )

    summary = [
        {
            "Sample": args.sample,
            "Variant_type": kind,
            "Evidence_flag": flag,
            "Observation_count": count,
            "Validation_status": "UNVALIDATED",
        }
        for (kind, flag), count in sorted(counts.items())
    ]

    output_dir = args.output_dir.expanduser().resolve()
    candidate_path = output_dir / f"{args.sample}_Preliminary_Candidates.csv"
    summary_path = output_dir / f"{args.sample}_Candidate_Evidence_Summary.csv"

    if candidate_path.exists() or summary_path.exists():
        raise FileExistsError("Refusing to overwrite existing candidate results")

    output_dir.mkdir(parents=True, exist_ok=True)

    for path, columns, rows in [
        (candidate_path, COLUMNS, candidates),
        (
            summary_path,
            [
                "Sample", "Variant_type", "Evidence_flag",
                "Observation_count", "Validation_status",
            ],
            summary,
        ),
    ]:
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            writer.writerows(rows)

    print("\nGENERALIZED SANGER CANDIDATE CONSOLIDATION")
    print("=" * 65)
    print("Sample:", args.sample)
    print("Raw SNV observations:", len(snvs))
    print("Distinct SNV alleles:", len(used_snvs))
    print("Raw deletion observations:", len(used_deletions))
    print("Raw insertion observations:", len(used_insertions))
    print("Total candidate observation rows:", len(candidates))

    print("\nEvidence summary:")
    for row in summary:
        print(
            f"  {row['Variant_type']} | "
            f"{row['Evidence_flag']}: {row['Observation_count']}"
        )

    print("\nSaved:")
    print(" ", candidate_path)
    print(" ", summary_path)
    print("\nAll observations remain UNVALIDATED.")


if __name__ == "__main__":
    main()
