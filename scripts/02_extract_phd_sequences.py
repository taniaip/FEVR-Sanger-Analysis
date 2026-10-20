#!/usr/bin/env python3

"""
Extract Sanger base calls and Phred scores from PHD files.

Inputs:
    AB1/PHD inventory CSV from 01_inventory_ab1_phd.py

Outputs:
    <sample>_All_Sanger_Reads.fastq
    <sample>_All_Sanger_Reads.fasta
    <sample>_KIF11_Genomic_Reads.fastq
    <sample>_MAP3K6_Plasmid_Reads.fastq
    <sample>_Initial_Read_QC.csv

No original sequencing files are modified.
"""

import argparse
import csv
import statistics
from pathlib import Path


def parse_phd(path):
    """Extract sequence and Phred quality from BEGIN_DNA section."""

    sequence = []
    qualities = []
    inside_dna = False

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()

            if line == "BEGIN_DNA":
                inside_dna = True
                continue

            if line == "END_DNA":
                break

            if not inside_dna or not line:
                continue

            fields = line.split()

            if len(fields) < 2:
                raise ValueError(
                    f"Malformed PHD base-call line in {path}: {line}"
                )

            base = fields[0].upper()

            if base not in {"A", "C", "G", "T", "N"}:
                raise ValueError(
                    f"Unexpected base {base!r} in {path}"
                )

            quality = int(fields[1])

            if quality < 0 or quality > 93:
                raise ValueError(
                    f"Invalid Phred score {quality} in {path}"
                )

            sequence.append(base)
            qualities.append(quality)

    if not sequence:
        raise ValueError(
            f"No PHD base calls found in {path}"
        )

    return "".join(sequence), qualities


def qc_metrics(sequence, qualities, min_length, min_mean_q, min_q20_fraction):
    length = len(sequence)
    mean_quality = statistics.mean(qualities)
    q20_fraction = sum(q >= 20 for q in qualities) / length
    n_fraction = sequence.count("N") / length

    if length < min_length:
        flag = "SHORT_READ"
    elif mean_quality < min_mean_q or q20_fraction < min_q20_fraction:
        flag = "LOW_MEAN_QUALITY"
    else:
        flag = "PASS_INITIAL_QC"

    return {
        "Read_length": length,
        "Mean_Phred_quality": round(mean_quality, 4),
        "Q20_fraction": round(q20_fraction, 4),
        "N_fraction": round(n_fraction, 4),
        "Initial_QC_flag": flag,
    }


def write_fastq(handle, read_id, sequence, qualities):
    encoded = "".join(chr(q + 33) for q in qualities)
    handle.write(f"@{read_id}\n{sequence}\n+\n{encoded}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Extract Sanger PHD sequences and initial QC."
    )

    parser.add_argument("--inventory", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)

    parser.add_argument("--min-length", type=int, default=100)
    parser.add_argument("--min-mean-q", type=float, default=20)
    parser.add_argument("--min-q20-fraction", type=float, default=0.5)

    args = parser.parse_args()

    inventory = args.inventory.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()

    if not inventory.is_file():
        parser.error(f"Inventory not found: {inventory}")

    if args.min_length < 1:
        parser.error("--min-length must be positive")

    if not 0 <= args.min_q20_fraction <= 1:
        parser.error("--min-q20-fraction must be between 0 and 1")

    with inventory.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        parser.error("Inventory is empty")

    selected = [
        row for row in rows
        if row["Sample"] == args.sample
    ]

    if not selected:
        parser.error(
            f"No reads for sample {args.sample}"
        )

    if len({row["Read_ID"] for row in selected}) != len(selected):
        parser.error("Duplicate Read_ID values detected")

    extracted = []
    qc_rows = []

    for row in selected:
        phd = Path(row["PHD_path"])

        if not phd.is_file():
            parser.error(f"Missing PHD file: {phd}")

        sequence, qualities = parse_phd(phd)

        metrics = qc_metrics(
            sequence,
            qualities,
            args.min_length,
            args.min_mean_q,
            args.min_q20_fraction,
        )

        extracted.append((row, sequence, qualities))

        qc_rows.append({
            "Sample": row["Sample"],
            "Read_ID": row["Read_ID"],
            "Sequencing_group": row["Sequencing_group"],
            "Gene": row["Gene"],
            "Target": row["Target"],
            "Direction": row["Direction"],
            "Filename": row["Filename"],
            **metrics,
        })

    output_dir.mkdir(parents=True, exist_ok=True)

    all_fastq = output_dir / f"{args.sample}_All_Sanger_Reads.fastq"
    all_fasta = output_dir / f"{args.sample}_All_Sanger_Reads.fasta"
    genomic_fastq = output_dir / f"{args.sample}_KIF11_Genomic_Reads.fastq"
    plasmid_fastq = output_dir / f"{args.sample}_MAP3K6_Plasmid_Reads.fastq"
    qc_csv = output_dir / f"{args.sample}_Initial_Read_QC.csv"

    with (
        all_fastq.open("w") as fq,
        all_fasta.open("w") as fa,
        genomic_fastq.open("w") as genomic,
        plasmid_fastq.open("w") as plasmid,
    ):
        for row, sequence, qualities in extracted:
            read_id = row["Read_ID"]

            write_fastq(fq, read_id, sequence, qualities)
            fa.write(f">{read_id}\n{sequence}\n")

            if row["Sequencing_group"] == "GENOMIC_FEVR":
                write_fastq(
                    genomic, read_id, sequence, qualities
                )

            elif row["Sequencing_group"] == "PLASMID":
                write_fastq(
                    plasmid, read_id, sequence, qualities
                )

    with qc_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(qc_rows[0].keys()),
        )
        writer.writeheader()
        writer.writerows(qc_rows)

    from collections import Counter

    counts = Counter(
        row["Initial_QC_flag"]
        for row in qc_rows
    )

    total_bases = sum(
        row["Read_length"]
        for row in qc_rows
    )

    print("\nSANGER PHD EXTRACTION — INITIAL QC")
    print("=" * 65)
    print("Sample:", args.sample)
    print("Reads extracted:", len(extracted))
    print("Total bases:", total_bases)

    print("\nInitial QC:")
    for flag in [
        "PASS_INITIAL_QC",
        "LOW_MEAN_QUALITY",
        "SHORT_READ",
    ]:
        print(f"  {flag}: {counts[flag]}")

    print("\nSaved:")
    for path in [
        all_fastq,
        all_fasta,
        genomic_fastq,
        plasmid_fastq,
        qc_csv,
    ]:
        print(" ", path)

    print("\nOriginal sequencing files were not modified.")


if __name__ == "__main__":
    main()
