#!/usr/bin/env python3
"""Run exploratory AB1 trace analysis for SNVs in an existing candidate CSV."""

import argparse
import csv
import subprocess
import sys
import tempfile
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--reverse-retain-bases", type=int)
    args = parser.parse_args()

    root = args.results_dir.expanduser().resolve()
    reference = args.reference.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()

    candidates_path = (
        root / "candidates" / f"{args.sample}_Preliminary_Candidates.csv"
    )

    with candidates_path.open(newline="", encoding="utf-8-sig") as handle:
        candidates = list(csv.DictReader(handle))

    variants = sorted({
        (
            row["Chromosome"],
            int(row["Position_GRCh38_1based"]),
            row["REF"].upper(),
            row["ALT"].upper(),
        )
        for row in candidates
        if row["Variant_type"] == "SNV"
        and len(row["REF"]) == 1
        and len(row["ALT"]) == 1
        and row["REF"].upper() in "ACGT"
        and row["ALT"].upper() in "ACGT"
        and row["REF"].upper() != row["ALT"].upper()
    })

    print("TRACE-AWARE CANDIDATE BATCH")
    print("=" * 70)
    print("Sample:", args.sample)
    print("Unique SNV candidates:", len(variants))
    print("Indels and non-SNV observations excluded.")

    if not variants:
        print("No SNV candidates to analyze.")
        return

    module = Path(__file__).with_name(
        "08c_generalized_trace_snv_evidence.py"
    )

    summary = []

    with tempfile.TemporaryDirectory() as temporary:
        temporary_dir = Path(temporary)

        for chromosome, position, ref, alt in variants:
            label = f"{chromosome}_{position}_{ref}_{alt}"
            temporary_output = temporary_dir / f"{label}.csv"

            command = [
                sys.executable,
                str(module),
                "--sample", args.sample,
                "--results-dir", str(root),
                "--reference", str(reference),
                "--chromosome", chromosome,
                "--position", str(position),
                "--ref", ref,
                "--alt", alt,
                "--output", str(temporary_output),
            ]

            if args.reverse_retain_bases is not None:
                command.extend([
                    "--reverse-retain-bases",
                    str(args.reverse_retain_bases),
                ])

            result = subprocess.run(
                command,
                text=True,
                capture_output=True,
            )

            status = "PASS" if result.returncode == 0 else "FAILED"

            summary.append({
                "Sample": args.sample,
                "Chromosome": chromosome,
                "Position_GRCh38_1based": position,
                "REF": ref,
                "ALT": alt,
                "Status": status,
                "Error": (
                    result.stderr.strip().splitlines()[-1]
                    if result.returncode and result.stderr.strip()
                    else ""
                ),
            })

            print(f"{label}: {status}")

            if result.returncode:
                print("  ", summary[-1]["Error"])
            else:
                output_dir.mkdir(parents=True, exist_ok=True)
                destination = output_dir / f"{args.sample}_{label}_Trace.csv"

                if destination.exists():
                    print("  Existing output preserved:", destination)
                else:
                    destination.write_bytes(temporary_output.read_bytes())

    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / f"{args.sample}_Trace_Batch_Summary.csv"

    if summary_path.exists():
        print("Existing summary preserved:", summary_path)
    else:
        with summary_path.open("x", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=list(summary[0]),
            )
            writer.writeheader()
            writer.writerows(summary)

    print("\nSummary:", summary_path)
    print("All observations remain UNVALIDATED.")


if __name__ == "__main__":
    main()
