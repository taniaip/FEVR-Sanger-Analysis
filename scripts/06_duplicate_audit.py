#!/usr/bin/env python3

"""
Audit Sanger AB1 files for byte-identical duplicates.

Uses SHA-256 file hashes and optionally joins genomic alignment QC.

A unique AB1 hash indicates a distinct file, not necessarily an
independent biological sample, library, or sequencing experiment.

Inputs:
    Inventory CSV from 01_inventory_ab1_phd.py
    Alignment QC CSV from 04_align_genomic_reads.py (optional)

Output:
    <sample>_AB1_Duplicate_Audit.csv
"""

import argparse
import csv
import hashlib
from collections import Counter, defaultdict
from pathlib import Path


def sha256_file(path):
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def load_csv(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def main():
    parser = argparse.ArgumentParser(
        description="Audit AB1 files for byte-identical duplicates."
    )

    parser.add_argument("--inventory", required=True, type=Path)
    parser.add_argument("--alignment-qc", type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)

    args = parser.parse_args()

    inventory_path = args.inventory.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()

    if not inventory_path.is_file():
        parser.error(f"Inventory not found: {inventory_path}")

    inventory = [
        row for row in load_csv(inventory_path)
        if row["Sample"] == args.sample
    ]

    if not inventory:
        parser.error(f"No inventory rows for sample {args.sample}")

    read_ids = [row["Read_ID"] for row in inventory]

    if len(read_ids) != len(set(read_ids)):
        parser.error("Duplicate Read_ID values in inventory")

    alignment_flags = {}

    if args.alignment_qc:
        qc_path = args.alignment_qc.expanduser().resolve()

        if not qc_path.is_file():
            parser.error(f"Alignment QC not found: {qc_path}")

        for row in load_csv(qc_path):
            if row["Sample"] != args.sample:
                continue

            read_id = row["Read_ID"]

            if read_id in alignment_flags:
                parser.error(
                    f"Duplicate alignment QC record: {read_id}"
                )

            alignment_flags[read_id] = row["Alignment_QC_flag"]

    hashes = {}
    hash_to_ids = defaultdict(list)

    for row in inventory:
        read_id = row["Read_ID"]
        ab1_path = Path(row["AB1_path"])

        if not ab1_path.is_file():
            parser.error(f"AB1 file missing: {ab1_path}")

        digest = sha256_file(ab1_path)

        hashes[read_id] = digest
        hash_to_ids[digest].append(read_id)

    output_rows = []

    for row in inventory:
        read_id = row["Read_ID"]
        digest = hashes[read_id]

        matching_ids = sorted(hash_to_ids[digest])
        other_ids = [
            identifier
            for identifier in matching_ids
            if identifier != read_id
        ]

        qc_flag = alignment_flags.get(read_id, "NOT_ASSESSED")

        # This is only a technical screening flag.
        # It does not establish biological independence.
        unique_passing_alignment = (
            qc_flag == "PASS"
            and len(matching_ids) == 1
        )

        output_rows.append({
            "Sample": row["Sample"],
            "Read_ID": read_id,
            "Sequencing_group": row["Sequencing_group"],
            "Gene": row["Gene"],
            "Target": row["Target"],
            "Direction": row["Direction"],
            "Filename": row["Filename"],
            "AB1_SHA256": digest,
            "Identical_AB1_count": len(matching_ids) - 1,
            "Duplicate_read_IDs": ";".join(other_ids),
            "Alignment_QC_flag": qc_flag,
            "Independent_passing_alignment_evidence": (
                unique_passing_alignment
            ),
        })

    output_dir.mkdir(parents=True, exist_ok=True)

    output_path = (
        output_dir / f"{args.sample}_AB1_Duplicate_Audit.csv"
    )

    with output_path.open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(output_rows[0].keys()),
        )
        writer.writeheader()
        writer.writerows(output_rows)

    duplicate_groups = [
        read_group
        for read_group in hash_to_ids.values()
        if len(read_group) > 1
    ]

    group_counts = Counter(
        row["Sequencing_group"]
        for row in inventory
    )

    passing_unique = sum(
        row["Independent_passing_alignment_evidence"]
        for row in output_rows
    )

    print("\nSANGER AB1 DUPLICATE AUDIT")
    print("=" * 65)
    print("Sample:", args.sample)
    print("AB1 files examined:", len(inventory))
    print("Distinct SHA-256 hashes:", len(hash_to_ids))
    print("Byte-identical duplicate groups:", len(duplicate_groups))
    print(
        "Files belonging to duplicate groups:",
        sum(len(group) for group in duplicate_groups)
    )

    print("\nSequencing groups:")
    for group, count in sorted(group_counts.items()):
        print(f"  {group}: {count}")

    print("\nUnique-file passing genomic alignments:", passing_unique)

    if duplicate_groups:
        print("\nDuplicate read groups:")
        for group in duplicate_groups:
            print(" ", ", ".join(sorted(group)))

    print("\nSaved:")
    print(" ", output_path)

    print(
        "\nNote: Distinct AB1 files do not establish "
        "independent biological evidence."
    )
    print("Original sequencing files were not modified.")


if __name__ == "__main__":
    main()
