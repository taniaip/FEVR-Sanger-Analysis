#!/usr/bin/env python3
"""Generate a generalized FEVR Sanger AB1/PHD inventory."""

import argparse
import csv
from collections import Counter
from pathlib import Path

from filename_parser import parse_filename
from metadata_overrides import load_overrides


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--overrides", type=Path)
    args = parser.parse_args()

    input_dir = args.input.expanduser().resolve()
    output = args.output.expanduser().resolve()

    if not input_dir.is_dir():
        parser.error(f"Input directory not found: {input_dir}")

    overrides = (
        load_overrides(args.overrides.expanduser())
        if args.overrides else {}
    )

    files = sorted(
        p for p in input_dir.rglob("*")
        if p.is_file()
        and p.suffix.lower() == ".ab1"
        and "Analysis" not in p.relative_to(input_dir).parts
    )

    if not files:
        parser.error("No AB1 files found")

    rows = []

    for index, ab1 in enumerate(files, start=1):
        phd = ab1.with_suffix(".phd.1")
        metadata = parse_filename(ab1.name, args.sample)

        key = (args.sample, ab1.name)

        if key in overrides:
            verified = overrides[key]
            metadata.update({
                "Gene": verified["Gene"],
                "Target": verified["Target"],
                "Direction": verified["Direction"],
                "Primer": verified["Direction"],
                "Parse_status": "VERIFIED_OVERRIDE",
            })

        rows.append({
            "Sample": args.sample,
            "Read_ID": f"read_{index:03d}",
            "Sequencing_group": metadata["Sequencing_group"],
            "Gene": metadata["Gene"],
            "Target": metadata["Target"],
            "Direction": metadata["Direction"],
            "Primer": metadata["Primer"],
            "Original_gene_label": metadata["Original_gene_label"],
            "Sample_preparation": metadata["Sample_preparation"],
            "Construct": metadata["Construct"],
            "Replicate": metadata["Replicate"],
            "Parse_status": metadata["Parse_status"],
            "Plate_well": "",
            "Run_ID": "",
            "Filename": ab1.name,
            "AB1_path": str(ab1.resolve()),
            "PHD_path": str(phd.resolve()),
            "PHD_exists": phd.is_file(),
            "Inventory_flag": (
                "MATCHED" if phd.is_file() else "MISSING_PHD"
            ),
        })

    output.parent.mkdir(parents=True, exist_ok=True)

    if output.exists():
        parser.error(f"Refusing to overwrite existing inventory: {output}")

    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    counts = Counter(row["Parse_status"] for row in rows)

    print("\nGENERALIZED FEVR INVENTORY")
    print("=" * 50)
    print("Sample:", args.sample)
    print("AB1 files:", len(rows))
    print("Matched PHD files:", sum(r["PHD_exists"] for r in rows))
    print("Parsing:", dict(counts))
    print("Saved:", output)
    print("Original sequencing files were not modified.")


if __name__ == "__main__":
    main()
