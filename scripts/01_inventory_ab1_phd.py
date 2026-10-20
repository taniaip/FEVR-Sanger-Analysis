#!/usr/bin/env python3

"""
Create an inventory of Sanger AB1 files and matching PHD files.

Supports genomic FEVR sequencing and MAP3K6 plasmid sequencing.

Example:
    python3 scripts/01_inventory_ab1_phd.py \
        --input ~/Desktop/FEVR-Project/426 \
        --sample 426 \
        --output inventory.csv

No original sequencing files are modified.
"""

import argparse
import csv
import re
from collections import Counter
from pathlib import Path


def parse_filename(filename, sample):
    """Extract available metadata from an AB1 filename."""

    stem = Path(filename).stem

    plasmid_match = re.search(
        rf"MAP3K6_+plasmid_+{re.escape(sample)}_+"
        r"(M13F|M13R|F[1-6]|R[1-6])(?=_|$)",
        stem,
        flags=re.IGNORECASE,
    )

    if plasmid_match:
        primer = plasmid_match.group(1).upper()

        return {
            "Sequencing_group": "PLASMID",
            "Gene": "MAP3K6",
            "Target": "PLASMID",
            "Direction": (
                "F" if primer.startswith("F") or primer == "M13F"
                else "R"
            ),
            "Primer": primer,
        }

    gene_match = re.search(
        r"(CDH5|KIF11|LRP5|NDP|TSPAN12)",
        stem,
        flags=re.IGNORECASE,
    )

    gene = gene_match.group(1).upper() if gene_match else ""

    exon_match = re.search(
        r"(?<![A-Za-z0-9])E(\d+)(?![A-Za-z0-9])",
        stem,
        flags=re.IGNORECASE,
    )

    target = f"E{exon_match.group(1)}" if exon_match else ""

    direction_match = re.search(
        r"(?:^|_)(F|R)(?:_|$)",
        stem,
        flags=re.IGNORECASE,
    )

    direction = (
        direction_match.group(1).upper()
        if direction_match else ""
    )

    return {
        "Sequencing_group": (
            "GENOMIC_FEVR" if gene else "UNCLASSIFIED"
        ),
        "Gene": gene,
        "Target": target,
        "Direction": direction,
        "Primer": direction,
    }


def build_inventory(input_dir, sample):
    """Collect AB1 files and check matching PHD files."""

    ab1_files = sorted(
        (
            path for path in input_dir.rglob("*")
            if path.is_file()
            and path.suffix.lower() == ".ab1"
            and "Analysis" not in path.relative_to(input_dir).parts
        ),
        key=lambda path: str(path.relative_to(input_dir)),
    )

    records = []

    for index, ab1 in enumerate(ab1_files, start=1):
        phd = ab1.with_suffix(".phd.1")
        metadata = parse_filename(ab1.name, sample)

        well_match = re.search(
            r"_([A-H]\d{2})_",
            ab1.stem,
            flags=re.IGNORECASE,
        )

        run_match = re.search(
            r"_(\d{3})$",
            ab1.stem,
        )

        records.append({
            "Sample": sample,
            "Read_ID": f"read_{index:03d}",
            **metadata,
            "Plate_well": (
                well_match.group(1).upper()
                if well_match else ""
            ),
            "Run_ID": (
                run_match.group(1)
                if run_match else ""
            ),
            "Filename": ab1.name,
            "AB1_path": str(ab1.resolve()),
            "PHD_path": str(phd.resolve()),
            "PHD_exists": phd.is_file(),
            "Inventory_flag": (
                "MATCHED"
                if phd.is_file()
                else "MISSING_PHD"
            ),
        })

    return records


def main():
    parser = argparse.ArgumentParser(
        description="Inventory Sanger AB1 and PHD files."
    )

    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Directory containing sequencing files.",
    )

    parser.add_argument(
        "--sample",
        required=True,
        help="Sample identifier.",
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output inventory CSV.",
    )

    args = parser.parse_args()

    input_dir = args.input.expanduser().resolve()
    output = args.output.expanduser().resolve()

    if not input_dir.is_dir():
        parser.error(f"Input directory not found: {input_dir}")

    records = build_inventory(input_dir, args.sample)

    if not records:
        parser.error("No AB1 files found.")

    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(records[0].keys()),
        )

        writer.writeheader()
        writer.writerows(records)

    counts = Counter(
        row["Sequencing_group"]
        for row in records
    )

    matched = sum(
        row["PHD_exists"]
        for row in records
    )

    print("\nSANGER AB1/PHD INVENTORY")
    print("=" * 60)
    print(f"Sample: {args.sample}")
    print(f"AB1 files: {len(records)}")
    print(f"Matched PHD files: {matched}")
    print(f"Missing PHD files: {len(records) - matched}")

    print("\nSequencing groups:")
    for group, count in sorted(counts.items()):
        print(f"  {group}: {count}")

    print(f"\nSaved: {output}")
    print("Original sequencing files were not modified.")


if __name__ == "__main__":
    main()
