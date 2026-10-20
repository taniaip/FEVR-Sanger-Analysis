#!/usr/bin/env python3

"""
Align trimmed genomic Sanger reads to a reference FASTA.

Requires:
    minimap2
    samtools

Outputs:
    <sample>_genomic.sam
    <sample>_genomic.sorted.bam
    <sample>_genomic.sorted.bam.bai
    <sample>_Read_Alignment_QC.csv

This script assesses read alignment, not validated variants.
Plasmid reads must not be aligned to a genomic reference.
"""

import argparse
import csv
import re
import shutil
import subprocess
from pathlib import Path


def run(command, stdout=None):
    subprocess.run(command, check=True, stdout=stdout)


def read_fastq_ids(path):
    ids = set()

    with path.open() as handle:
        while True:
            header = handle.readline()

            if not header:
                break

            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()

            if not sequence or not plus or not quality:
                raise ValueError("Incomplete FASTQ record")

            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError("Invalid FASTQ record")

            read_id = header[1:].strip().split()[0]

            if read_id in ids:
                raise ValueError(f"Duplicate FASTQ ID: {read_id}")

            ids.add(read_id)

    return ids


def cigar_metrics(cigar):
    operations = re.findall(r"(\d+)([MIDNSHP=X])", cigar)

    aligned_query = sum(
        int(length)
        for length, op in operations
        if op in "MI=X"
    )

    reference_length = sum(
        int(length)
        for length, op in operations
        if op in "MDN=X"
    )

    return aligned_query, reference_length


def parse_sam(sam_path):
    alignments = {}

    with sam_path.open() as handle:
        for line in handle:
            if line.startswith("@"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 11:
                raise ValueError("Malformed SAM record")

            read_id = fields[0]
            flag = int(fields[1])

            if flag & 0x100 or flag & 0x800:
                continue

            if read_id in alignments:
                raise ValueError(
                    f"Multiple primary alignments: {read_id}"
                )

            chromosome = fields[2]
            start = int(fields[3])
            mapq = int(fields[4])
            cigar = fields[5]

            unmapped = bool(flag & 0x4)

            if unmapped:
                alignments[read_id] = {
                    "Chromosome": "",
                    "Start_GRCh38_1based": "",
                    "End_GRCh38_1based": "",
                    "Alignment_strand": "",
                    "MAPQ": 0,
                    "CIGAR": "",
                    "Aligned_query_bases": 0,
                }
                continue

            aligned_query, reference_length = cigar_metrics(cigar)

            alignments[read_id] = {
                "Chromosome": chromosome,
                "Start_GRCh38_1based": start,
                "End_GRCh38_1based": start + reference_length - 1,
                "Alignment_strand": "-" if flag & 16 else "+",
                "MAPQ": mapq,
                "CIGAR": cigar,
                "Aligned_query_bases": aligned_query,
            }

    return alignments


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--fastq", required=True, type=Path)
    parser.add_argument("--inventory", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument(
        "--reference-build",
        required=True,
        choices=["GRCh38"],
        help="Declared reference build. Only GRCh38 is supported.",
    )
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--min-aligned-fraction", type=float, default=0.8)

    args = parser.parse_args()

    for program in ["minimap2", "samtools"]:
        if shutil.which(program) is None:
            parser.error(f"{program} is not installed or not in PATH")

    fastq = args.fastq.expanduser().resolve()
    inventory_path = args.inventory.expanduser().resolve()
    reference = args.reference.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()

    for path in [fastq, inventory_path, reference]:
        if not path.is_file():
            parser.error(f"File not found: {path}")

    if not 0 <= args.min_aligned_fraction <= 1:
        parser.error("Invalid minimum aligned fraction")

    with inventory_path.open(newline="") as handle:
        inventory = list(csv.DictReader(handle))

    genomic_metadata = {
        row["Read_ID"]: row
        for row in inventory
        if row["Sample"] == args.sample
        and row["Sequencing_group"] == "GENOMIC_FEVR"
    }

    read_ids = read_fastq_ids(fastq)

    if not read_ids:
        parser.error("No reads found in FASTQ")

    if not read_ids.issubset(genomic_metadata):
        parser.error(
            "FASTQ contains reads not identified as genomic FEVR reads"
        )

    output_dir.mkdir(parents=True, exist_ok=True)

    sam_path = output_dir / f"{args.sample}_genomic.sam"
    bam_path = output_dir / f"{args.sample}_genomic.sorted.bam"
    qc_path = output_dir / f"{args.sample}_Read_Alignment_QC.csv"

    command = [
        "minimap2",
        "-a",
        "-x", "map-ont",
        "--secondary=no",
        str(reference),
        str(fastq),
    ]

    print("Running minimap2...", flush=True)

    with sam_path.open("w") as handle:
        run(command, stdout=handle)

    print("Sorting BAM...", flush=True)

    run([
        "samtools", "sort",
        "-o", str(bam_path),
        str(sam_path),
    ])

    run(["samtools", "index", str(bam_path)])

    alignments = parse_sam(sam_path)

    rows = []

    for read_id in sorted(read_ids):
        metadata = genomic_metadata[read_id]
        alignment = alignments.get(read_id)

        if alignment is None:
            raise RuntimeError(
                f"No SAM alignment record for {read_id}"
            )

        # Count bases directly from the trimmed FASTQ.
        # The SAM SEQ field may be reverse-complemented.
        with fastq.open() as handle:
            for line in handle:
                if line.strip() == f"@{read_id}":
                    sequence = handle.readline().strip()
                    break
            else:
                raise RuntimeError(f"Read not found: {read_id}")

        trimmed_length = len(sequence)
        aligned_query = alignment["Aligned_query_bases"]

        aligned_fraction = (
            aligned_query / trimmed_length
            if trimmed_length else 0
        )

        if not alignment["Chromosome"]:
            qc_flag = "UNMAPPED"
        elif alignment["MAPQ"] < args.min_mapq:
            qc_flag = "LOW_MAPQ"
        elif aligned_fraction < args.min_aligned_fraction:
            qc_flag = "LOW_ALIGNED_FRACTION"
        else:
            qc_flag = "PASS"

        rows.append({
            "Sample": args.sample,
            "Read_ID": read_id,
            "Sequencing_group": metadata["Sequencing_group"],
            "Gene": metadata["Gene"],
            "Target": metadata["Target"],
            "Direction": metadata["Direction"],
            "Trimmed_length": trimmed_length,
            **alignment,
            "Aligned_fraction": round(aligned_fraction, 4),
            "Alignment_QC_flag": qc_flag,
        })

    with qc_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
        )
        writer.writeheader()
        writer.writerows(rows)

    from collections import Counter

    counts = Counter(row["Alignment_QC_flag"] for row in rows)

    print("\nGENOMIC SANGER ALIGNMENT QC")
    print("=" * 65)
    print("Sample:", args.sample)
    print("Genomic reads aligned:", len(rows))

    for flag, count in sorted(counts.items()):
        print(f"  {flag}: {count}")

    print("\nRead-level alignment:")
    for row in rows:
        print(
            f"  {row['Read_ID']} "
            f"{row['Chromosome']}:"
            f"{row['Start_GRCh38_1based']}-"
            f"{row['End_GRCh38_1based']} "
            f"MAPQ={row['MAPQ']} "
            f"fraction={row['Aligned_fraction']} "
            f"{row['Alignment_QC_flag']}"
        )

    print("\nSaved:")
    print(" ", sam_path)
    print(" ", bam_path)
    print(" ", str(bam_path) + ".bai")
    print(" ", qc_path)


if __name__ == "__main__":
    main()
