#!/usr/bin/env python3
"""Audit unvalidated Sanger insertion observations across matching amplicon reads."""

import argparse
import csv
import re
import subprocess
from collections import Counter
from pathlib import Path

COLUMNS = [
    "Sample", "Gene", "Target", "Read_ID", "Direction",
    "Chromosome", "Position_GRCh38_1based", "REF", "ALT",
    "Variant_type", "Insertion_length",
    "Left_flank_quality", "Right_flank_quality", "Both_flanks_Q20",
    "Exact_insertion_support_reads", "Other_insertion_reads",
    "No_insertion_Q20_reads", "No_insertion_lowQ_reads",
    "Not_fully_covered_reads", "Evidence_flag", "Validation_status",
]


def read_csv(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def cigar_ops(cigar):
    matches = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    if not matches or "".join(f"{n}{op}" for n, op in matches) != cigar:
        raise ValueError(f"Invalid CIGAR: {cigar}")
    return [(int(n), op) for n, op in matches]


def parse_alignment(fields):
    read_id = fields[0]
    flag = int(fields[1])
    chrom = fields[2]
    position = int(fields[3])
    sequence = fields[9].upper()
    quality = fields[10]

    if flag & (4 | 256 | 2048):
        raise ValueError(f"Unexpected SAM flag: {read_id}")
    if sequence == "*" or quality == "*" or len(sequence) != len(quality):
        raise ValueError(f"Missing/mismatched sequence: {read_id}")

    scores = [ord(char) - 33 for char in quality]
    bases = {}
    insertions = []
    q = 0
    r = position

    for length, op in cigar_ops(fields[5]):
        if op in {"M", "=", "X"}:
            for offset in range(length):
                bases[r + offset] = scores[q + offset]
            q += length
            r += length

        elif op == "I":
            insertions.append((r - 1, sequence[q:q + length]))
            q += length

        elif op in {"D", "N"}:
            r += length

        elif op == "S":
            q += length

        elif op in {"H", "P"}:
            pass

    if q != len(sequence):
        raise ValueError(f"CIGAR query length mismatch: {read_id}")

    return {
        "chrom": chrom,
        "bases": bases,
        "insertions": insertions,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
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
        if row["Sample"] == args.sample and row["Alignment_QC_flag"] == "PASS":
            read_id = row["Read_ID"]
            if read_id in passing:
                raise ValueError(f"Duplicate QC read: {read_id}")
            passing[read_id] = row

    if not passing:
        parser.error("No passing reads found")

    raw = [
        row for row in read_csv(raw_path)
        if row["Sample"] == args.sample and row["Variant_type"] == "INSERTION"
    ]

    if not raw:
        output_dir = args.output_dir.expanduser().resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"{args.sample}_Insertion_Evidence_Audit.csv"
        if output_path.exists():
            raise FileExistsError(f"Refusing to overwrite: {output_path}")
        with output_path.open("w", newline="", encoding="utf-8") as handle:
            csv.DictWriter(handle, fieldnames=COLUMNS).writeheader()
        print("No insertion observations; empty evidence CSV saved.")
        print("Saved:", output_path)
        return

    for row in raw:
        if row["Read_ID"] not in passing:
            raise ValueError(f"Insertion read failed QC: {row['Read_ID']}")

    sam = subprocess.run(
        ["samtools", "view", "-F", "2308", str(bam)],
        capture_output=True, text=True, check=True,
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
        raise RuntimeError("Passing BAM reads and QC records do not match")

    results = []

    for candidate in raw:
        chrom = candidate["Chromosome"]
        anchor = int(candidate["Position_GRCh38_1based"])
        ref = candidate["REF"].upper()
        alt = candidate["ALT"].upper()

        if len(ref) != 1 or not alt.startswith(ref) or len(alt) <= len(ref):
            raise ValueError(f"Unexpected insertion representation: {ref}>{alt}")

        inserted = alt[len(ref):]

        if len(inserted) != int(candidate["Indel_length"]):
            raise ValueError("Insertion length mismatch")

        left = candidate["Left_flank_quality"]
        right = candidate["Right_flank_quality"]
        left_q = int(left) if left else None
        right_q = int(right) if right else None

        both_q20 = (
            left_q is not None and right_q is not None
            and left_q >= args.min_base_q
            and right_q >= args.min_base_q
        )

        exact = []
        other = []
        intact_q20 = []
        intact_lowq = []
        not_covered = []

        for read_id, alignment in sorted(alignments.items()):
            metadata = passing[read_id]

            if (
                metadata["Gene"] != candidate["Gene"]
                or metadata["Target"] != candidate["Target"]
            ):
                continue

            if alignment["chrom"] != chrom:
                not_covered.append(read_id)
                continue

            site_insertions = [
                sequence
                for position, sequence in alignment["insertions"]
                if position == anchor
            ]

            if inserted in site_insertions:
                exact.append(read_id)
                continue

            if site_insertions:
                other.append(read_id)
                continue

            bases = alignment["bases"]

            if anchor not in bases or anchor + 1 not in bases:
                not_covered.append(read_id)
                continue

            if (
                bases[anchor] >= args.min_base_q
                and bases[anchor + 1] >= args.min_base_q
            ):
                intact_q20.append(read_id)
            else:
                intact_lowq.append(read_id)

        if not both_q20:
            evidence_flag = "LOW_FLANK_QUALITY"
        elif len(exact) < 2:
            evidence_flag = "Q20_INSERTION_UNCONFIRMED"
        else:
            evidence_flag = "MULTIREAD_INSERTION_REQUIRES_VALIDATION"

        results.append({
            "Sample": args.sample,
            "Gene": candidate["Gene"],
            "Target": candidate["Target"],
            "Read_ID": candidate["Read_ID"],
            "Direction": candidate["Direction"],
            "Chromosome": chrom,
            "Position_GRCh38_1based": anchor,
            "REF": ref,
            "ALT": alt,
            "Variant_type": "INSERTION",
            "Insertion_length": len(inserted),
            "Left_flank_quality": left,
            "Right_flank_quality": right,
            "Both_flanks_Q20": both_q20,
            "Exact_insertion_support_reads": ";".join(exact),
            "Other_insertion_reads": ";".join(other),
            "No_insertion_Q20_reads": ";".join(intact_q20),
            "No_insertion_lowQ_reads": ";".join(intact_lowq),
            "Not_fully_covered_reads": ";".join(not_covered),
            "Evidence_flag": evidence_flag,
            "Validation_status": "UNVALIDATED",
        })

    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / f"{args.sample}_Insertion_Evidence_Audit.csv"

    if output.exists():
        raise FileExistsError(f"Refusing to overwrite: {output}")

    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(results)

    print("\nGENERALIZED SANGER INSERTION EVIDENCE")
    print("=" * 65)
    print("Sample:", args.sample)
    print("Passing reads examined:", len(alignments))
    print("Insertion observations:", len(results))
    print("Minimum flank quality:", args.min_base_q)

    print("\nEvidence categories:")
    for flag, count in sorted(Counter(r["Evidence_flag"] for r in results).items()):
        print(f"  {flag}: {count}")

    print("\nSaved:", output)
    print("All insertions remain UNVALIDATED and unnormalized.")
    print("Exact read support does not establish biological independence.")


if __name__ == "__main__":
    main()
