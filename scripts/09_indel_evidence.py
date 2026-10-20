#!/usr/bin/env python3

"""
Assess preliminary Sanger deletion observations.

For each raw deletion:
- Check Phred quality of the two flanking bases.
- Count reads with the exact same anchored deletion.
- Count other deletions overlapping the deleted interval.
- Count reads with an intact reference-spanning alignment,
  separating Q20 and lower-quality spanning evidence.
- Identify reads that do not fully cover the deleted interval.

Deletion coordinates and representations are not normalized.
Evidence is technical read-level evidence, not a genotype or
independent biological confirmation.

Requires samtools.
"""

import argparse
import csv
import re
import subprocess
from collections import Counter
from pathlib import Path


COLUMNS = [
    "Sample", "Gene", "Target", "Read_ID", "Direction",
    "Chromosome", "Position_GRCh38_1based", "REF", "ALT",
    "Variant_type", "Deletion_length", "Left_flank_quality",
    "Right_flank_quality", "Both_flanks_Q20",
    "Exact_deletion_support_reads",
    "Overlapping_deletion_reads",
    "Intact_reference_Q20_reads",
    "Intact_reference_lowQ_reads",
    "Not_fully_covered_reads",
    "Evidence_flag", "Validation_status",
]


def read_csv(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def cigar_ops(cigar):
    ops = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    if not ops or "".join(f"{n}{op}" for n, op in ops) != cigar:
        raise ValueError(f"Invalid CIGAR: {cigar}")
    return [(int(n), op) for n, op in ops]


def parse_alignment(fields):
    read_id = fields[0]
    flag = int(fields[1])
    chrom = fields[2]
    pos = int(fields[3])
    seq = fields[9]
    qual = fields[10]

    if flag & (4 | 256 | 2048):
        raise ValueError(f"Unexpected alignment flag: {read_id}")
    if seq == "*" or qual == "*" or len(seq) != len(qual):
        raise ValueError(f"Missing/mismatched sequence: {read_id}")

    scores = [ord(c) - 33 for c in qual]
    bases = {}
    deletions = []
    q = 0
    r = pos

    for length, op in cigar_ops(fields[5]):
        if op in {"M", "=", "X"}:
            for offset in range(length):
                bases[r + offset] = scores[q + offset]
            r += length
            q += length

        elif op == "D":
            deletions.append((r, r + length - 1))
            r += length

        elif op == "N":
            r += length

        elif op in {"I", "S"}:
            q += length

        elif op in {"H", "P"}:
            pass

    if q != len(seq):
        raise ValueError(f"CIGAR query length mismatch: {read_id}")

    return {
        "chrom": chrom,
        "bases": bases,
        "deletions": deletions,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Audit Sanger deletion flank quality and read evidence."
    )
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--alignment-qc", required=True, type=Path)
    parser.add_argument("--raw-indels", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--min-base-q", type=int, default=20)
    args = parser.parse_args()

    if not 0 <= args.min_base_q <= 93:
        parser.error("--min-base-q must be between 0 and 93")

    paths = [
        args.bam.expanduser().resolve(),
        args.alignment_qc.expanduser().resolve(),
        args.raw_indels.expanduser().resolve(),
    ]

    for path in paths:
        if not path.is_file():
            parser.error(f"Missing input: {path}")

    bam, qc_path, raw_path = paths

    passing = {}
    for row in read_csv(qc_path):
        if (
            row["Sample"] == args.sample
            and row["Alignment_QC_flag"] == "PASS"
        ):
            read_id = row["Read_ID"]
            if read_id in passing:
                raise ValueError(f"Duplicate QC read: {read_id}")
            passing[read_id] = row

    if not passing:
        parser.error("No passing reads found")

    raw = [
        row for row in read_csv(raw_path)
        if row["Sample"] == args.sample
    ]

    if not raw:
        parser.error("No raw indel observations found")

    for row in raw:
        if row["Variant_type"] != "DELETION":
            parser.error(
                "This audit currently supports deletions only; "
                "insertion evidence requires a separate assessment"
            )
        if row["Read_ID"] not in passing:
            parser.error(
                f"Indel read did not pass alignment QC: {row['Read_ID']}"
            )

    sam = subprocess.run(
        ["samtools", "view", "-F", "2308", str(bam)],
        capture_output=True,
        text=True,
        check=True,
    )

    alignments = {}
    for line in sam.stdout.splitlines():
        fields = line.split("\t")
        read_id = fields[0]
        if read_id not in passing:
            continue
        if read_id in alignments:
            raise ValueError(f"Multiple primary alignments: {read_id}")
        alignments[read_id] = parse_alignment(fields)

    if set(alignments) != set(passing):
        raise RuntimeError("Passing BAM reads and QC table do not match")

    output_rows = []

    for candidate in raw:
        chrom = candidate["Chromosome"]
        anchor_pos = int(candidate["Position_GRCh38_1based"])
        ref = candidate["REF"].upper()
        alt = candidate["ALT"].upper()

        if not ref.startswith(alt) or len(ref) <= len(alt):
            raise ValueError(
                f"Unexpected deletion representation: {ref}>{alt}"
            )

        deletion_length = len(ref) - len(alt)
        deletion_start = anchor_pos + len(alt)
        deletion_end = deletion_start + deletion_length - 1

        if deletion_length != int(candidate["Indel_length"]):
            raise ValueError("Deletion length mismatch")

        left = candidate["Left_flank_quality"]
        right = candidate["Right_flank_quality"]

        left_q = int(left) if left not in ("", None) else None
        right_q = int(right) if right not in ("", None) else None

        both_q20 = (
            left_q is not None
            and right_q is not None
            and left_q >= args.min_base_q
            and right_q >= args.min_base_q
        )

        exact = []
        overlapping = []
        intact_q20 = []
        intact_lowq = []
        not_covered = []

        for read_id, alignment in sorted(alignments.items()):
            if alignment["chrom"] != chrom:
                not_covered.append(read_id)
                continue

            deletions = alignment["deletions"]

            if (deletion_start, deletion_end) in deletions:
                exact.append(read_id)
                continue

            if any(
                d_start <= deletion_end and d_end >= deletion_start
                for d_start, d_end in deletions
            ):
                overlapping.append(read_id)
                continue

            # Require the anchor, every deleted reference position,
            # and the immediately following reference base.
            required = range(anchor_pos, deletion_end + 2)
            bases = alignment["bases"]

            if not all(pos in bases for pos in required):
                not_covered.append(read_id)
                continue

            if all(bases[pos] >= args.min_base_q for pos in required):
                intact_q20.append(read_id)
            else:
                intact_lowq.append(read_id)

        # A single observation is never considered independently
        # confirmed merely because its flanking bases are Q20.
        if not both_q20:
            evidence_flag = "LOW_FLANK_QUALITY"
        elif len(exact) < 2:
            evidence_flag = "Q20_DELETION_UNCONFIRMED"
        else:
            evidence_flag = "MULTIREAD_DELETION_REQUIRES_VALIDATION"

        output_rows.append({
            "Sample": args.sample,
            "Gene": candidate["Gene"],
            "Target": candidate["Target"],
            "Read_ID": candidate["Read_ID"],
            "Direction": candidate["Direction"],
            "Chromosome": chrom,
            "Position_GRCh38_1based": anchor_pos,
            "REF": ref,
            "ALT": alt,
            "Variant_type": "DELETION",
            "Deletion_length": deletion_length,
            "Left_flank_quality": left if left != "" else "",
            "Right_flank_quality": right if right != "" else "",
            "Both_flanks_Q20": both_q20,
            "Exact_deletion_support_reads": ";".join(exact),
            "Overlapping_deletion_reads": ";".join(overlapping),
            "Intact_reference_Q20_reads": ";".join(intact_q20),
            "Intact_reference_lowQ_reads": ";".join(intact_lowq),
            "Not_fully_covered_reads": ";".join(not_covered),
            "Evidence_flag": evidence_flag,
            "Validation_status": "UNVALIDATED",
        })

    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / f"{args.sample}_Indel_Quality_Evidence_Audit.csv"

    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(output_rows)

    counts = Counter(row["Evidence_flag"] for row in output_rows)

    print("\nSANGER INDEL QUALITY AND EVIDENCE AUDIT")
    print("=" * 65)
    print("Sample:", args.sample)
    print("Passing reads examined:", len(alignments))
    print("Raw deletion observations:", len(output_rows))
    print("Minimum flank quality:", args.min_base_q)

    print("\nEvidence categories:")
    for flag, count in sorted(counts.items()):
        print(f"  {flag}: {count}")

    print("\nQ20-flanked deletion observations:")
    found = False
    for row in output_rows:
        if row["Both_flanks_Q20"]:
            found = True
            print(
                f"  {row['Read_ID']} "
                f"{row['Chromosome']}:"
                f"{row['Position_GRCh38_1based']} "
                f"{row['REF']}>{row['ALT']} "
                f"flanks={row['Left_flank_quality']}/"
                f"{row['Right_flank_quality']} "
                f"{row['Evidence_flag']}"
            )

    if not found:
        print("  None")

    print("\nSaved:")
    print(" ", output)
    print("\nAll deletion observations remain UNVALIDATED.")
    print("Raw indel representations have not been normalized.")


if __name__ == "__main__":
    main()
