#!/usr/bin/env python3
"""Run the generalized FEVR pipeline sequentially across eligible samples."""

import csv
import subprocess
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
INPUT_ROOT = Path.home() / "Desktop" / "FEVR-Project"
REFERENCE = INPUT_ROOT / "references" / "FEVR_GRCh38_5chrom.fa"
RESULTS = ROOT / "results"
LOGS = RESULTS / "batch_logs"
SUMMARY = LOGS / "batch_status.csv"

COMPLETED_SAMPLES = {"428", "476", "551"}

def file_counts(folder):
    ab1 = [
        p for p in folder.rglob("*")
        if p.is_file() and p.suffix.lower() == ".ab1"
    ]
    phd = [
        p for p in folder.rglob("*")
        if p.is_file() and p.name.lower().endswith(".phd.1")
    ]

    ab1_names = Counter(p.stem.casefold() for p in ab1)
    phd_names = Counter(p.name[:-6].casefold() for p in phd)
    missing = sum((ab1_names - phd_names).values())

    return len(ab1), missing

def main():
    if not REFERENCE.is_file():
        raise SystemExit(f"Missing reference: {REFERENCE}")

    if not Path(str(REFERENCE) + ".fai").is_file():
        raise SystemExit("Missing reference FASTA index.")

    if not (ROOT / "run_generalized_pipeline.py").is_file():
        raise SystemExit("Missing generalized pipeline runner.")

    samples = []

    for folder in sorted(INPUT_ROOT.iterdir()):
        if not folder.is_dir() or folder.name == "references":
            continue

        count, missing = file_counts(folder)

        if count:
            samples.append((folder, count, missing))

    print("FEVR GENERALIZED BATCH PLAN")
    print("=" * 75)

    for folder, count, missing in samples:
        output = RESULTS / f"{folder.name}_generalized_batch"

        if folder.name in COMPLETED_SAMPLES:
            status = "PREVIOUSLY_ANALYZED"
        elif missing:
            status = "BLOCKED_MISSING_PHD"
        elif output.exists():
            status = "OUTPUT_EXISTS"
        else:
            status = "READY"

        print(f"{folder.name:<10} {status:<24} AB1={count}")

    if "--run" not in sys.argv:
        print("\nDRY RUN ONLY: No analyses were started.")
        print("Use --run to execute the batch.")
        return

    LOGS.mkdir(parents=True, exist_ok=True)

    rows = []

    for folder, count, missing in samples:
        sample = folder.name
        output = RESULTS / f"{sample}_generalized_batch"

        if sample in COMPLETED_SAMPLES:
            status = "PREVIOUSLY_ANALYZED"
            returncode = ""
        elif missing:
            status = "BLOCKED_MISSING_PHD"
            returncode = ""
        elif output.exists():
            status = "OUTPUT_EXISTS"
            returncode = ""
        else:
            command = [
                sys.executable,
                str(ROOT / "run_generalized_pipeline.py"),
                "--sample", sample,
                "--input", str(folder),
                "--reference", str(REFERENCE),
                "--reference-build", "GRCh38",
                "--output-dir", str(output),
            ]

            log_path = LOGS / f"{sample}.log"

            print(f"\nRunning sample {sample}...", flush=True)

            with log_path.open("w", encoding="utf-8") as log:
                log.write("Command: " + " ".join(command) + "\n")
                log.flush()

                try:
                    result = subprocess.run(
                        command,
                        cwd=ROOT,
                        stdout=log,
                        stderr=subprocess.STDOUT,
                        check=False,
                    )
                    returncode = result.returncode
                    status = "COMPLETED" if returncode == 0 else "FAILED"
                except Exception as exc:
                    status = "ERROR"
                    returncode = ""
                    log.write(f"\nBatch execution error: {exc}\n")

            print(f"Sample {sample}: {status}", flush=True)

        rows.append({
            "Sample": sample,
            "AB1_count": count,
            "Missing_PHD": missing,
            "Status": status,
            "Return_code": returncode,
            "Output_directory": str(output),
            "Timestamp": datetime.now().isoformat(timespec="seconds"),
        })

        with SUMMARY.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)

    print("\nBatch finished.")
    print("Status summary:", SUMMARY)
    print("Logs:", LOGS)

if __name__ == "__main__":
    main()
