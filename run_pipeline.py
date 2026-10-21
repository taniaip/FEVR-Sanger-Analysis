#!/usr/bin/env python3
"""Run the existing preliminary KIF11 E10/MAP3K6 Sanger scripts in order."""
import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
NAMES = ["01_inventory_ab1_phd.py", "02_extract_phd_sequences.py",
         "03_quality_trim.py", "04_align_genomic_reads.py",
         "05_alignment_coverage.py", "06_duplicate_audit.py",
         "07_extract_variants.py", "08_snv_evidence.py",
         "09_indel_evidence.py", "10_consolidate_candidates.py",
         "11_plasmid_primer_qc.py", "12_generate_workbook.py"]

def main():
    p = argparse.ArgumentParser(description=__doc__)
    for flag in ("sample", "input", "reference", "output-dir"):
        p.add_argument("--" + flag, required=True)
    a = p.parse_args()
    s = a.sample
    source = Path(a.input).expanduser().resolve()
    ref = Path(a.reference).expanduser().resolve()
    out = Path(a.output_dir).expanduser().resolve()
    if not s.isdigit() or not source.is_dir() or not ref.is_file():
        p.error("Check numeric sample ID, input directory, and reference FASTA.")
    if out.exists():
        p.error(f"Refusing to overwrite existing directory: {out}")
    if not list(source.glob("*.ab1")) and not list(source.glob("*.AB1")):
        p.error("No AB1 files found.")
    for name in NAMES:
        if not (ROOT / "scripts" / name).is_file():
            p.error(f"Missing script: {name}")
    for name in ("minimap2", "samtools"):
        if not shutil.which(name):
            p.error(f"Missing dependency: {name}")
    out.mkdir(parents=True)
    def path(folder, name):
        return out / folder / (s + "_" + name)
    inventory = out / (s + "_Inventory.csv")
    bam = path("alignment", "genomic.sorted.bam")
    qc = path("alignment", "Read_Alignment_QC.csv")
    def run(n, *args):
        print(f"\\nSTEP {n:02d}/12: {NAMES[n-1]}", flush=True)
        subprocess.run([sys.executable, str(ROOT / "scripts" / NAMES[n-1]),
                        *map(str, args)], check=True, cwd=ROOT)
    try:
        run(1, "--input", source, "--sample", s, "--output", inventory)
        with inventory.open(newline="", encoding="utf-8-sig") as h:
            rows = list(csv.DictReader(h))
        if not rows or any(r.get("Sample") != s for r in rows):
            raise ValueError("Empty inventory or incorrect sample.")
        genomic = [r for r in rows if r.get("Sequencing_group") == "GENOMIC_FEVR"]
        plasmid = [r for r in rows if r.get("Sequencing_group") == "PLASMID"]
        if not genomic or not plasmid or len(genomic) + len(plasmid) != len(rows):
            raise ValueError("Unsupported sequencing groups or missing read group.")
        if any(r.get("Gene") != "KIF11" or r.get("Target") != "E10" for r in genomic):
            raise ValueError("Only KIF11 E10 genomic reads are supported.")
        if any(r.get("Gene") != "MAP3K6" for r in plasmid):
            raise ValueError("Only MAP3K6 plasmid reads are supported.")
        if any("MISSING" in str(r).upper() for r in rows):
            raise ValueError("Inventory indicates missing sequencing files.")
        print(f"Validated {len(genomic)} KIF11 E10 and {len(plasmid)} MAP3K6 reads.", flush=True)
        run(2, "--inventory", inventory, "--sample", s, "--output-dir", out/"initial")
        run(3, "--fastq", path("initial", "All_Sanger_Reads.fastq"),
            "--inventory", inventory, "--sample", s, "--output-dir", out/"trimmed")
        run(4, "--fastq", path("trimmed", "KIF11_trimmed.fastq"),
            "--inventory", inventory, "--reference", ref, "--reference-build", "GRCh38",
            "--sample", s, "--output-dir", out/"alignment")
        run(5, "--bam", bam, "--alignment-qc", qc, "--sample", s,
            "--gene", "KIF11", "--target", "E10", "--output-dir", out/"coverage")
        run(6, "--inventory", inventory, "--alignment-qc", qc,
            "--sample", s, "--output-dir", out/"duplicates")
        run(7, "--bam", bam, "--alignment-qc", qc, "--reference", ref,
            "--sample", s, "--output-dir", out/"variants")
        run(8, "--bam", bam, "--alignment-qc", qc,
            "--raw-snvs", path("variants", "Preliminary_SNV_Read_Evidence.csv"),
            "--sample", s, "--output-dir", out/"variants")
        run(9, "--bam", bam, "--alignment-qc", qc,
            "--raw-indels", path("variants", "Preliminary_Indel_Read_Evidence.csv"),
            "--sample", s, "--output-dir", out/"variants")
        run(10, "--raw-snvs", path("variants", "Preliminary_SNV_Read_Evidence.csv"),
            "--snv-evidence", path("variants", "SNV_Evidence_Audit.csv"),
            "--raw-indels", path("variants", "Preliminary_Indel_Read_Evidence.csv"),
            "--indel-evidence", path("variants", "Indel_Quality_Evidence_Audit.csv"),
            "--sample", s, "--output-dir", out/"candidates")
        run(11, "--inventory", inventory,
            "--trim-qc", path("trimmed", "Quality_Trim_Summary.csv"),
            "--sample", s, "--output-dir", out/"plasmid")
        workbook = out / (s + "_FEVR_Sanger_Analysis.xlsx")
        run(12, "--results-dir", out, "--sample", s, "--gene", "KIF11",
            "--target", "E10", "--reference-label", "GRCh38 chr10 (unverified)",
            "--output", workbook)
        if not workbook.is_file():
            raise RuntimeError("Final workbook missing.")
        print(f"\\nCOMPLETED: {workbook}\\nPreliminary, unvalidated results.")
    except (subprocess.CalledProcessError, ValueError, RuntimeError) as e:
        print(f"PIPELINE STOPPED: {e}\\nPartial outputs: {out}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
