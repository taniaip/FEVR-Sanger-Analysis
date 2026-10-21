#!/usr/bin/env python3
"""Multi-gene Sanger alignment coverage; not exon completeness or variant validation."""

import argparse
import csv
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

DETAIL_FIELDS = [
    "Sample", "Read_ID", "Gene", "Target", "Chromosome",
    "Primer_direction", "Direction_parser_accepted", "Alignment_strand",
    "Start_GRCh38_1based", "End_GRCh38_1based", "CIGAR",
    "Aligned_reference_bases", "Insertion_bases", "Deletion_bases",
    "Softclipped_bases", "Alignment_QC_flag", "Review_flags",
]

COVERAGE_FIELDS = [
    "Sample", "Gene", "Target", "Chromosome",
    "Position_GRCh38_1based", "Read_depth",
    "Supporting_read_IDs", "Parser_forward_reads",
    "Parser_reverse_reads", "Unknown_direction_reads",
    "Both_parser_directions",
]

SUMMARY_FIELDS = [
    "Sample", "Gene", "Target", "Total_trimmed_reads",
    "Passing_reads", "Excluded_reads", "Unmapped_reads",
    "Passing_parser_forward_reads", "Passing_parser_reverse_reads",
    "Passing_unknown_direction_reads", "Covered_reference_positions",
    "Positions_with_both_parser_directions", "Maximum_read_depth",
    "Mean_depth_at_covered_positions",
    "Exon_completeness_assessed",
]


def parse_cigar(cigar):
    operations = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    if not operations or "".join(n + op for n, op in operations) != cigar:
        raise ValueError(f"Invalid CIGAR: {cigar}")
    return [(int(n), op) for n, op in operations]


def load_qc(path, sample):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {
            "Sample", "Read_ID", "Gene", "Target", "Direction",
            "Parse_status", "Chromosome", "Alignment_QC_flag",
        }
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Alignment QC is missing required columns")

        rows = {}
        for row in reader:
            if row["Sample"] != sample:
                continue
            read_id = row["Read_ID"]
            if read_id in rows:
                raise ValueError(f"Duplicate QC Read_ID: {read_id}")
            rows[read_id] = row

    if not rows:
        raise ValueError(f"No QC records for sample {sample}")
    return rows


def direction_class(row):
    if row["Parse_status"] == "OK" and row["Direction"] in {"F", "R"}:
        return row["Direction"]
    return "UNKNOWN"


def process_alignment(fields, metadata, coverage):
    read_id = fields[0]
    flag = int(fields[1])
    chromosome = fields[2]
    reference_position = int(fields[3])
    cigar = fields[5]

    if flag & (0x4 | 0x100 | 0x800):
        raise ValueError(f"Unexpected SAM flag for {read_id}")
    if chromosome != metadata["Chromosome"]:
        raise ValueError(f"Chromosome mismatch for {read_id}")
    if cigar != metadata["CIGAR"]:
        raise ValueError(f"CIGAR mismatch for {read_id}")

    start = reference_position
    aligned_bases = 0
    insertion_bases = 0
    deletion_bases = 0
    softclipped_bases = 0

    for length, operation in parse_cigar(cigar):
        if operation in {"M", "=", "X"}:
            for position in range(reference_position, reference_position + length):
                coverage[(metadata["Gene"], metadata["Target"], chromosome, position)].add(read_id)
            aligned_bases += length
            reference_position += length
        elif operation == "I":
            insertion_bases += length
        elif operation == "D":
            deletion_bases += length
            reference_position += length
        elif operation == "N":
            reference_position += length
        elif operation == "S":
            softclipped_bases += length

    direction = direction_class(metadata)
    return {
        "Sample": metadata["Sample"],
        "Read_ID": read_id,
        "Gene": metadata["Gene"],
        "Target": metadata["Target"],
        "Chromosome": chromosome,
        "Primer_direction": metadata["Direction"],
        "Direction_parser_accepted": direction != "UNKNOWN",
        "Alignment_strand": "-" if flag & 16 else "+",
        "Start_GRCh38_1based": start,
        "End_GRCh38_1based": reference_position - 1,
        "CIGAR": cigar,
        "Aligned_reference_bases": aligned_bases,
        "Insertion_bases": insertion_bases,
        "Deletion_bases": deletion_bases,
        "Softclipped_bases": softclipped_bases,
        "Alignment_QC_flag": metadata["Alignment_QC_flag"],
        "Review_flags": metadata.get("Review_flags", ""),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--alignment-qc", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    if not shutil.which("samtools"):
        parser.error("samtools not found in PATH")

    bam = args.bam.expanduser().resolve()
    qc_path = args.alignment_qc.expanduser().resolve()
    output = args.output_dir.expanduser().resolve()

    for path in (bam, qc_path):
        if not path.is_file():
            parser.error(f"Missing input: {path}")

    qc = load_qc(qc_path, args.sample)
    passing = {
        read_id: row for read_id, row in qc.items()
        if row["Alignment_QC_flag"] == "PASS"
    }

    prefix = args.sample
    details_path = output / f"{prefix}_Generalized_Alignment_Details.csv"
    coverage_path = output / f"{prefix}_Generalized_Position_Coverage.csv"
    summary_path = output / f"{prefix}_Generalized_Target_Summary.csv"

    for path in (details_path, coverage_path, summary_path):
        if path.exists():
            parser.error(f"Output exists; refusing overwrite: {path}")

    result = subprocess.run(
        ["samtools", "view", "-F", "2308", str(bam)],
        capture_output=True, text=True, check=True,
    )

    coverage = defaultdict(set)
    details = []
    seen = set()

    for line in result.stdout.splitlines():
        fields = line.split("\t")
        read_id = fields[0]

        if read_id not in passing:
            continue
        if read_id in seen:
            raise ValueError(f"Multiple primary alignments: {read_id}")

        seen.add(read_id)
        details.append(process_alignment(fields, passing[read_id], coverage))

    missing = set(passing) - seen
    if missing:
        raise RuntimeError(f"Passing reads missing from BAM: {sorted(missing)}")

    coverage_rows = []
    positions_by_target = defaultdict(list)

    for (gene, target, chromosome, position), read_ids in sorted(coverage.items()):
        directions = Counter(direction_class(qc[read_id]) for read_id in read_ids)
        f_count = directions["F"]
        r_count = directions["R"]
        unknown_count = directions["UNKNOWN"]

        row = {
            "Sample": args.sample,
            "Gene": gene,
            "Target": target,
            "Chromosome": chromosome,
            "Position_GRCh38_1based": position,
            "Read_depth": len(read_ids),
            "Supporting_read_IDs": ";".join(sorted(read_ids)),
            "Parser_forward_reads": f_count,
            "Parser_reverse_reads": r_count,
            "Unknown_direction_reads": unknown_count,
            "Both_parser_directions": f_count > 0 and r_count > 0,
        }
        coverage_rows.append(row)
        positions_by_target[(gene, target)].append(row)

    target_groups = defaultdict(list)
    for row in qc.values():
        target_groups[(row["Gene"], row["Target"])].append(row)

    summary_rows = []
    for (gene, target), rows in sorted(target_groups.items()):
        passed = [row for row in rows if row["Alignment_QC_flag"] == "PASS"]
        positions = positions_by_target[(gene, target)]
        direction_counts = Counter(direction_class(row) for row in passed)
        depths = [row["Read_depth"] for row in positions]

        summary_rows.append({
            "Sample": args.sample,
            "Gene": gene,
            "Target": target,
            "Total_trimmed_reads": len(rows),
            "Passing_reads": len(passed),
            "Excluded_reads": len(rows) - len(passed),
            "Unmapped_reads": sum(row["Alignment_QC_flag"] == "UNMAPPED" for row in rows),
            "Passing_parser_forward_reads": direction_counts["F"],
            "Passing_parser_reverse_reads": direction_counts["R"],
            "Passing_unknown_direction_reads": direction_counts["UNKNOWN"],
            "Covered_reference_positions": len(positions),
            "Positions_with_both_parser_directions": sum(
                row["Both_parser_directions"] for row in positions
            ),
            "Maximum_read_depth": max(depths, default=0),
            "Mean_depth_at_covered_positions": round(
                sum(depths) / len(depths), 3
            ) if depths else 0,
            "Exon_completeness_assessed": False,
        })

    output.mkdir(parents=True, exist_ok=True)

    for path, fields, rows in (
        (details_path, DETAIL_FIELDS, details),
        (coverage_path, COVERAGE_FIELDS, coverage_rows),
        (summary_path, SUMMARY_FIELDS, summary_rows),
    ):
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)

    print("\nGENERALIZED SANGER COVERAGE")
    print("=" * 65)
    print("Sample:", args.sample)
    print("Total trimmed reads:", len(qc))
    print("Passing reads included:", len(details))
    print("Genes:", len({row["Gene"] for row in summary_rows}))
    print("Gene/target groups:", len(summary_rows))
    print("Reference-position records:", len(coverage_rows))
    print("\nSaved:")
    for path in (details_path, coverage_path, summary_path):
        print(" ", path)
    print("\nExon completeness and variant validation were not assessed.")


if __name__ == "__main__":
    main()
