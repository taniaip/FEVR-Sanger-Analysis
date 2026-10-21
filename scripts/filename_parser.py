#!/usr/bin/env python3
"""Parse FEVR Sanger filenames without modifying original files."""

import re
from pathlib import Path

GENE_ALIASES = {"KIFF11": "KIF11"}

KNOWN_GENES = {
    "ADAM10", "CDH5", "FZD4", "FZD7", "GAPDH",
    "IQGAP1", "KIF11", "LRP5", "MAP3K6", "NDP",
    "OLFM3", "TSPAN12", "ZNF408",
}


def normalize_gene(label):
    label = label.upper()
    return GENE_ALIASES.get(label, label)


def parse_filename(filename, sample):
    stem = Path(filename).stem
    parts = stem.split("_")
    original_gene = parts[0].upper()
    gene = normalize_gene(original_gene)
    sample = str(sample)

    result = {
        "Gene": gene if gene in KNOWN_GENES else "",
        "Original_gene_label": original_gene,
        "Target": "",
        "Direction": "",
        "Sample_preparation": "UNSPECIFIED",
        "Construct": "",
        "Sequencing_group": "UNCLASSIFIED",
        "Primer": "",
        "Replicate": "",
        "Parse_status": "REVIEW",
    }

    if gene not in KNOWN_GENES:
        return result

    if gene == "MAP3K6":
        result["Sequencing_group"] = "PLASMID"

        if re.search(rf"MUT2-{re.escape(sample)}(?:_|$)", stem, re.I):
            result["Construct"] = "MUT2"
        elif "plasmid" in stem.lower():
            result["Construct"] = "ORIGINAL"

        primer_match = re.search(
            r"(?:^|_)(M13F|M13R|F[1-6]|R[1-6])(?=_|$)",
            stem,
            re.I,
        )
        if primer_match:
            primer = primer_match.group(1).upper()
            result["Primer"] = primer
            result["Direction"] = (
                "F" if primer.startswith("F") or primer == "M13F"
                else "R"
            )

        result["Target"] = "PLASMID"
        result["Parse_status"] = (
            "OK" if result["Construct"] and result["Primer"]
            else "REVIEW"
        )
        return result

    result["Sequencing_group"] = "GENOMIC_FEVR"

    if re.search(r"WGA\.", stem, re.I):
        result["Sample_preparation"] = "WGA"
    elif re.search(r"cDNA", stem, re.I):
        result["Sample_preparation"] = "CDNA"
    elif re.search(r"Blood", stem, re.I):
        result["Sample_preparation"] = "BLOOD"

    sample_pattern = re.escape(sample)
    sample_match = re.search(
        rf"(?<![A-Za-z0-9])(?:cDNA|WGA)?{sample_pattern}(?:-(\d+))?(?![A-Za-z0-9])",
        stem,
        re.I,
    )
    if sample_match and sample_match.group(1):
        result["Replicate"] = sample_match.group(1)

    direction_match = re.search(
        r"(?:^|_)(F[1-9]?|R[1-9]?)(?=_|$)", stem, re.I
    )
    if direction_match:
        result["Primer"] = direction_match.group(1).upper()
        result["Direction"] = result["Primer"][0]

    # Recognize exon labels wherever they occur in the filename.
    exon_matches = [
        p.upper() for p in parts
        if re.fullmatch(
            r"E?\d+(?:E\d+|[A-Z]|-\d+)*",
            p,
            re.I,
        )
    ]

    # Exclude numeric sample identifiers and sequencing run identifiers.
    exon_matches = [
        p for p in exon_matches
        if p != sample.upper()
        and not re.fullmatch(r"\d{3}", p)
    ]

    if exon_matches:
        target = exon_matches[0]
        if re.fullmatch(r"\d+[A-Z]?", target):
            target = "E" + target

        target = re.sub(
            r"^E(\d+)",
            lambda match: "E" + match.group(1).zfill(2),
            target,
        )

        result["Target"] = target

    if gene == "NDP" and len(parts) > 1:
        if parts[1].upper() in {"A", "B", "C"}:
            result["Target"] = parts[1].upper()

    if gene == "GAPDH":
        result["Target"] = "CONTROL"

    # Some filenames lack sufficient information for automatic assignment.
    if result["Target"] and result["Direction"] and sample_match:
        result["Parse_status"] = "OK"

    return result
