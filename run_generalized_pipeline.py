#!/usr/bin/env python3
"""Run the 11-stage generalized FEVR Sanger pipeline through sample-level QC reporting."""

import argparse
import csv
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent

SCRIPTS = [
    "01b_generalized_inventory.py",
    "02_extract_phd_sequences.py",
    "03_quality_trim.py",
    "04b_generalized_alignment.py",
    "05b_generalized_coverage.py",
    "07b_generalized_variant_extraction.py",
    "08b_generalized_snv_evidence.py",
    "09b_generalized_deletion_evidence.py",
    "09c_generalized_insertion_evidence.py",
    "10b_generalized_candidates.py",
    "11b_generalized_qc_report.py",
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--reference-build", default="GRCh38")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--overrides", type=Path)
    args = parser.parse_args()

    sample = args.sample
    source = args.input.expanduser().resolve()
    reference = args.reference.expanduser().resolve()
    output = args.output_dir.expanduser().resolve()

    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", sample):
        parser.error("Sample ID must contain only letters, numbers, _ or -.")

    if not source.is_dir():
        parser.error(f"Input directory does not exist: {source}")

    if not reference.is_file():
        parser.error(f"Reference FASTA does not exist: {reference}")

    if output.exists():
        parser.error(f"Refusing to overwrite existing directory: {output}")

    if not any(
        path.is_file() and path.suffix.lower() == ".ab1"
        for path in source.iterdir()
    ):
        parser.error("No AB1 files found in the input directory.")

    if args.overrides:
        overrides = args.overrides.expanduser().resolve()
        if not overrides.is_file():
            parser.error(f"Overrides file does not exist: {overrides}")
    else:
        overrides = None

    for script in SCRIPTS:
        if not (ROOT / "scripts" / script).is_file():
            parser.error(f"Missing script: {script}")

    for dependency in ("minimap2", "samtools"):
        if not shutil.which(dependency):
            parser.error(f"Missing dependency: {dependency}")

    inventory = output / f"{sample}_Inventory.csv"
    initial_fastq = output / "initial" / f"{sample}_All_Sanger_Reads.fastq"
    genomic_fastq = output / "trimmed" / f"{sample}_KIF11_trimmed.fastq"

    alignment_dir = output / "alignment"
    bam = alignment_dir / f"{sample}_generalized_genomic.sorted.bam"
    alignment_qc = alignment_dir / f"{sample}_Generalized_Alignment_QC.csv"
    coverage_dir = output / "coverage"
    variants_dir = output / "variants"
    snv_dir = output / "snv_evidence"
    deletion_dir = output / "deletion_evidence"
    insertion_dir = output / "insertion_evidence"
    candidates_dir = output / "candidates"

    raw_snvs = variants_dir / f"{sample}_Preliminary_SNV_Read_Evidence.csv"
    raw_indels = variants_dir / f"{sample}_Preliminary_Indel_Read_Evidence.csv"
    snv_evidence = snv_dir / f"{sample}_SNV_Evidence_Audit.csv"
    deletion_evidence = deletion_dir / f"{sample}_Indel_Quality_Evidence_Audit.csv"
    insertion_evidence = insertion_dir / f"{sample}_Insertion_Evidence_Audit.csv"
    candidates = candidates_dir / f"{sample}_Preliminary_Candidates.csv"

    def run(stage, *arguments):
        script = SCRIPTS[stage - 1]
        print(
            f"\nSTEP {stage}/{len(SCRIPTS)}: {script}",
            flush=True,
        )
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "scripts" / script),
                *map(str, arguments),
            ],
            check=True,
            cwd=ROOT,
        )

    output.mkdir(parents=True, exist_ok=False)

    try:
        inventory_args = [
            "--input", source,
            "--sample", sample,
            "--output", inventory,
        ]

        if overrides:
            inventory_args.extend(["--overrides", overrides])

        run(1, *inventory_args)

        with inventory.open(newline="", encoding="utf-8-sig") as handle:
            rows = list(csv.DictReader(handle))

        if not rows:
            raise ValueError("Inventory is empty.")

        if any(row.get("Sample") != sample for row in rows):
            raise ValueError("Inventory contains an unexpected sample ID.")

        read_ids = [row["Read_ID"] for row in rows]

        if len(read_ids) != len(set(read_ids)):
            raise ValueError("Duplicate Read_ID values in inventory.")

        groups = {row.get("Sequencing_group") for row in rows}

        if not groups.issubset({"GENOMIC_FEVR", "PLASMID"}):
            raise ValueError(f"Unsupported sequencing groups: {groups}")

        genomic_count = sum(
            row["Sequencing_group"] == "GENOMIC_FEVR"
            for row in rows
        )

        if genomic_count == 0:
            raise ValueError("No genomic FEVR reads found.")

        print(
            f"Inventory: {len(rows)} reads; "
            f"{genomic_count} genomic reads.",
            flush=True,
        )

        run(
            2,
            "--inventory", inventory,
            "--sample", sample,
            "--output-dir", output / "initial",
        )

        run(
            3,
            "--fastq", initial_fastq,
            "--inventory", inventory,
            "--sample", sample,
            "--output-dir", output / "trimmed",
        )

        if not genomic_fastq.is_file() or genomic_fastq.stat().st_size == 0:
            raise ValueError("No retained genomic FASTQ reads.")

        run(
            4,
            "--fastq", genomic_fastq,
            "--inventory", inventory,
            "--reference", reference,
            "--reference-build", args.reference_build,
            "--sample", sample,
            "--output-dir", alignment_dir,
        )

        with alignment_qc.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames or not {
                "Sample", "Read_ID", "Alignment_QC_flag"
            }.issubset(reader.fieldnames):
                raise ValueError("Alignment QC is missing required columns.")
            qc_rows = list(reader)

        if not qc_rows:
            raise RuntimeError(
                "INSUFFICIENT SEQUENCING EVIDENCE: "
                "alignment QC contains no reads."
            )

        if any(row["Sample"] != sample for row in qc_rows):
            raise ValueError("Alignment QC contains an unexpected sample ID.")

        qc_read_ids = [row["Read_ID"] for row in qc_rows]
        if len(qc_read_ids) != len(set(qc_read_ids)):
            raise ValueError("Duplicate Read_ID values in alignment QC.")

        passing_count = sum(
            row["Alignment_QC_flag"] == "PASS"
            for row in qc_rows
        )

        print(
            f"Alignment QC: {passing_count}/{len(qc_rows)} "
            "genomic reads passed.",
            flush=True,
        )

        if passing_count == 0:
            raise RuntimeError(
                "INSUFFICIENT SEQUENCING EVIDENCE: "
                "zero genomic reads passed alignment QC. "
                "No variant-negative conclusion can be drawn."
            )

        run(
            5,
            "--bam", bam,
            "--alignment-qc", alignment_qc,
            "--sample", sample,
            "--output-dir", coverage_dir,
        )

        run(
            6,
            "--bam", bam,
            "--alignment-qc", alignment_qc,
            "--reference", reference,
            "--sample", sample,
            "--output-dir", variants_dir,
        )

        run(
            7,
            "--bam", bam,
            "--alignment-qc", alignment_qc,
            "--raw-snvs", raw_snvs,
            "--sample", sample,
            "--output-dir", snv_dir,
        )

        run(
            8,
            "--bam", bam,
            "--alignment-qc", alignment_qc,
            "--raw-indels", raw_indels,
            "--sample", sample,
            "--output-dir", deletion_dir,
        )

        run(
            9,
            "--bam", bam,
            "--alignment-qc", alignment_qc,
            "--raw-indels", raw_indels,
            "--sample", sample,
            "--output-dir", insertion_dir,
        )

        run(
            10,
            "--raw-snvs", raw_snvs,
            "--snv-evidence", snv_evidence,
            "--raw-indels", raw_indels,
            "--deletion-evidence", deletion_evidence,
            "--insertion-evidence", insertion_evidence,
            "--sample", sample,
            "--output-dir", candidates_dir,
        )

        def count_rows(csv_path):
            with csv_path.open(newline="", encoding="utf-8-sig") as handle:
                return sum(1 for _ in csv.DictReader(handle))

        expected = count_rows(raw_snvs) + count_rows(raw_indels)
        observed = count_rows(candidates)

        if expected != observed:
            raise RuntimeError(
                f"Candidate count mismatch: {expected} raw vs {observed} consolidated"
            )

        run(
            11,
            "--sample", sample,
            "--results-dir", output,
            "--output-dir", output / "qc_report",
        )

        print("\nGENERALIZED PIPELINE COMPLETED")
        print("=" * 65)
        print("Sample:", sample)
        print("Output:", output)
        print("Stages completed: inventory through sample-level QC reporting")
        print("QC reports:", output / "qc_report")
        print("Preliminary observation rows:", observed)
        print("All variant observations remain UNVALIDATED.")
        print("No genotype or clinical interpretation performed.")

    except (subprocess.CalledProcessError, ValueError, RuntimeError) as error:
        print(
            f"\nPIPELINE STOPPED: {error}\n"
            f"Partial outputs retained at: {output}",
            file=sys.stderr,
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
