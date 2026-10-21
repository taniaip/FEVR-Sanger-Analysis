#!/usr/bin/env python3
"""Load manually verified Sanger filename metadata overrides."""

import csv
from pathlib import Path

from filename_parser import KNOWN_GENES, normalize_gene

REQUIRED_COLUMNS = {
    "Sample",
    "Filename",
    "Verified_gene",
    "Verified_target",
    "Verified_direction",
    "Verification_source",
    "Review_notes",
}


def load_overrides(path):
    """Return verified overrides keyed by (sample, filename)."""
    path = Path(path)

    if not path.is_file():
        raise FileNotFoundError(f"Override file not found: {path}")

    overrides = {}
    seen = set()

    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)

        if not reader.fieldnames or not REQUIRED_COLUMNS.issubset(reader.fieldnames):
            raise ValueError("Override CSV is missing required columns")

        for line_number, row in enumerate(reader, start=2):
            sample = row["Sample"].strip()
            filename = row["Filename"].strip()
            key = (sample, filename)

            if not sample or not filename:
                raise ValueError(f"Missing sample or filename at line {line_number}")

            if key in seen:
                raise ValueError(f"Duplicate override record: {key}")
            seen.add(key)

            gene = row["Verified_gene"].strip().upper()
            target = row["Verified_target"].strip().upper()
            direction = row["Verified_direction"].strip().upper()
            source = row["Verification_source"].strip()

            # Unverified records remain pending.
            if not any((gene, target, direction, source)):
                continue

            if not all((gene, target, direction, source)):
                raise ValueError(
                    f"Incomplete verification at line {line_number}"
                )

            gene = normalize_gene(gene)

            if gene not in KNOWN_GENES:
                raise ValueError(f"Unknown verified gene at line {line_number}")

            if direction not in {"F", "R"}:
                raise ValueError(
                    f"Invalid verified direction at line {line_number}"
                )

            overrides[key] = {
                "Gene": gene,
                "Target": target,
                "Direction": direction,
                "Verification_source": source,
                "Review_notes": row["Review_notes"].strip(),
            }

    return overrides
