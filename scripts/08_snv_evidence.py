#!/usr/bin/env python3

"""
Audit preliminary Sanger SNVs against all passing genomic reads.

For each candidate allele, classify read observations as:
- Q20 alternate-base support
- Q20 reference-base support
- Q20 other-base support
- Low-quality base observation
- Deletion or reference-skip overlap

Only primary alignments passing the supplied alignment QC are used.

This is a read-evidence audit, not variant calling or genotyping.
"""

import argparse
import csv
import re
import subprocess
from collections import defaultdict
from pathlib import Path


COLUMNS = [
    "Sample", "Gene", "Target", "Chromosome",
    "Position_GRCh38_1based", "REF", "ALT", "Variant_type",
    "Independent_ALT_reads", "Independent_REF_reads",
    "Independent_other_base_reads", "Supporting_read_IDs",
    "Opposing_read_IDs", "Other_base_read_IDs",
    "Low_quality_read_IDs", "Deletion_or_skip_read_IDs",
    "ALT_primer_directions", "ALT_alignment_strands",
    "Evidence_flag", "Validation_status",
]


def load_csv(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def cigar_operations(cigar):
    operations = re.findall(r"(\d+)([MIDNSHP=X])", cigar)

    if not operations:
        raise ValueError(f"Invalid CIGAR: {cigar}")

    reconstructed = "".join(
        f"{length}{operation}"
        for length, operation in operations
    )

    if reconstructed != cigar:
        raise ValueError(f"Malformed CIGAR: {cigar}")

    return [(int(length), operation) for length, operation in operations]


def parse_alignment(fields):
    read_id = fields[0]
    flag = int(fields[1])
    chromosome = fields[2]
    start = int(fields[3])
    cigar = fields[5]
    sequence = fields[9].upper()
    quality_string = fields[10]

    if flag & (0x4 | 0x100 | 0x800):
        raise ValueError(f"Unexpected SAM flag: {read_id}")

    if sequence == "*" or quality_string == "*":
        raise ValueError(f"Missing sequence/quality: {read_id}")

    qualities = [ord(char) - 33 for char in quality_string]

    if len(sequence) != len(qualities):
        raise ValueError(f"Sequence/quality mismatch: {read_id}")

    observations = {}
    query_index = 0
    reference_position = start

    for length, operation in cigar_operations(cigar):

        if operation in {"M", "=", "X"}:
            for offset in range(length):
                observations[reference_position + offset] = (
                    sequence[query_index + offset],
                    qualities[query_index + offset],
                )

            query_index += length
            reference_position += length

        elif operation in {"D", "N"}:
            for position in range(
                reference_position,
                reference_position + length,
            ):
                observations[position] = ("DELETION_OR_SKIP", None)

            reference_position += length

        elif operation in {"I", "S"}:
            query_index += length

        elif operation in {"H", "P"}:
            pass

    if query_index != len(sequence):
        raise ValueError(f"CIGAR query mismatch: {read_id}")

    return {
        "Read_ID": read_id,
        "Chromosome": chromosome,
        "Strand": "-" if flag & 16 else "+",
        "Observations": observations,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Audit preliminary SNV evidence across passing reads."
    )

    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--alignment-qc", required=True, type=Path)
    parser.add_argument("--raw-snvs", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--min-base-q", type=int, default=20)

    args = parser.parse_args()

    if not 0 <= args.min_base_q <= 93:
        parser.error("--min-base-q must be between 0 and 93")

    bam = args.bam.expanduser().resolve()
    qc_path = args.alignment_qc.expanduser().resolve()
    snv_path = args.raw_snvs.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()

    for path in [bam, qc_path, snv_path]:
        if not path.is_file():
            parser.error(f"Missing input: {path}")

    passing = {}

    for row in load_csv(qc_path):
        if (
            row["Sample"] == args.sample
            and row["Alignment_QC_flag"] == "PASS"
        ):
            read_id = row["Read_ID"]

            if read_id in passing:
                raise ValueError(f"Duplicate QC read: {read_id}")

            passing[read_id] = row

    if not passing:
        parser.error("No passing reads in alignment QC")

    raw_snvs = [
        row for row in load_csv(snv_path)
        if row["Sample"] == args.sample
    ]

    if not raw_snvs:
        parser.error("No preliminary SNV observations")

    candidates = {}

    for row in raw_snvs:
        key = (
            row["Gene"],
            row["Target"],
            row["Chromosome"],
            int(row["Position_GRCh38_1based"]),
            row["REF"].upper(),
            row["ALT"].upper(),
        )
        candidates[key] = row

    result = subprocess.run(
        ["samtools", "view", "-F", "2308", str(bam)],
        capture_output=True,
        text=True,
        check=True,
    )

    alignments = {}

    for line in result.stdout.splitlines():
        fields = line.split("\t")
        read_id = fields[0]

        if read_id not in passing:
            continue

        if read_id in alignments:
            raise ValueError(
                f"Multiple primary alignments: {read_id}"
            )

        alignments[read_id] = parse_alignment(fields)

    missing = set(passing) - set(alignments)

    if missing:
        raise RuntimeError(
            f"Passing reads absent from BAM: {sorted(missing)}"
        )

    output_rows = []

    for key in sorted(candidates):
        gene, target, chromosome, position, ref, alt = key

        alt_ids = []
        ref_ids = []
        other_ids = []
        low_quality_ids = []
        deletion_ids = []

        alt_directions = set()
        alt_strands = set()

        for read_id, alignment in sorted(alignments.items()):
            if alignment["Chromosome"] != chromosome:
                continue

            observation = alignment["Observations"].get(position)

            if observation is None:
                continue

            base, quality = observation

            if base == "DELETION_OR_SKIP":
                deletion_ids.append(read_id)
                continue

            if quality < args.min_base_q or base not in "ACGT":
                low_quality_ids.append(read_id)
                continue

            if base == alt:
                alt_ids.append(read_id)
                alt_directions.add(passing[read_id]["Direction"])
                alt_strands.add(alignment["Strand"])

            elif base == ref:
                ref_ids.append(read_id)

            else:
                other_ids.append(read_id)

        if not alt_ids:
            evidence_flag = "NO_QUALIFYING_ALT_SUPPORT"

        elif ref_ids or other_ids:
            evidence_flag = "CONFLICTING_REF_ALT_EVIDENCE"

        else:
            evidence_flag = "ALT_SUPPORT_REQUIRES_VALIDATION"

        output_rows.append({
            "Sample": args.sample,
            "Gene": gene,
            "Target": target,
            "Chromosome": chromosome,
            "Position_GRCh38_1based": position,
            "REF": ref,
            "ALT": alt,
            "Variant_type": "SNV",
            "Independent_ALT_reads": len(alt_ids),
            "Independent_REF_reads": len(ref_ids),
            "Independent_other_base_reads": len(other_ids),
            "Supporting_read_IDs": ";".join(alt_ids),
            "Opposing_read_IDs": ";".join(ref_ids),
            "Other_base_read_IDs": ";".join(other_ids),
            "Low_quality_read_IDs": ";".join(low_quality_ids),
            "Deletion_or_skip_read_IDs": ";".join(deletion_ids),
            "ALT_primer_directions": ";".join(sorted(alt_directions)),
            "ALT_alignment_strands": ";".join(sorted(alt_strands)),
            "Evidence_flag": evidence_flag,
            "Validation_status": "UNVALIDATED",
        })

    output_dir.mkdir(parents=True, exist_ok=True)

    output_path = (
        output_dir / f"{args.sample}_SNV_Evidence_Audit.csv"
    )

    with output_path.open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(output_rows)

    from collections import Counter

    counts = Counter(
        row["Evidence_flag"] for row in output_rows
    )

    print("\nSANGER SNV EVIDENCE AUDIT")
    print("=" * 65)
    print("Sample:", args.sample)
    print("Passing reads examined:", len(alignments))
    print("Distinct SNV alleles:", len(output_rows))
    print("Minimum base quality:", args.min_base_q)

    print("\nEvidence categories:")

    for flag, count in sorted(counts.items()):
        print(f"  {flag}: {count}")

    print("\nCandidates with qualifying ALT support:")

    supported = [
        row for row in output_rows
        if row["Independent_ALT_reads"] > 0
    ]

    if not supported:
        print("  None")
    else:
        for row in supported:
            print(
                f"  {row['Chromosome']}:"
                f"{row['Position_GRCh38_1based']} "
                f"{row['REF']}>{row['ALT']} "
                f"ALT={row['Independent_ALT_reads']} "
                f"REF={row['Independent_REF_reads']} "
                f"OTHER={row['Independent_other_base_reads']} "
                f"{row['Evidence_flag']}"
            )

    print("\nSaved:")
    print(" ", output_path)

    print("\nAll candidates remain UNVALIDATED.")
    print("Read-level support does not establish biological independence.")


if __name__ == "__main__":
    main()
