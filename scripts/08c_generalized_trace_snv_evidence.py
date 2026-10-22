#!/usr/bin/env python3
"""Exploratory AB1 chromatogram evidence at specified genomic SNV sites.

This module does not call genotypes, validate variants, or override alignment QC.
Diagnostic trimming is optional and explicitly recorded.
"""

import argparse
import csv
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pysam
from Bio import SeqIO


COMPLEMENT = str.maketrans("ACGT", "TGCA")


def read_csv(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def align_sequence(sequence, read_id, reference):
    result = subprocess.run(
        ["minimap2", "-a", "-x", "sr", str(reference), "-"],
        input=f">{read_id}\n{sequence}\n",
        text=True,
        capture_output=True,
        check=True,
    )

    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "alignment.sam"
        path.write_text(result.stdout)

        with pysam.AlignmentFile(str(path), "r") as handle:
            alignments = [
                alignment
                for alignment in handle
                if not alignment.is_unmapped
                and not alignment.is_secondary
                and not alignment.is_supplementary
            ]

            if len(alignments) != 1:
                raise ValueError(
                    f"{read_id}: expected one primary mapped alignment; "
                    f"found {len(alignments)}"
                )

            alignment = alignments[0]

            return {
                "chromosome": alignment.reference_name,
                "reverse": alignment.is_reverse,
                "mapq": alignment.mapping_quality,
                "cigar": alignment.cigarstring,
                "aligned_pairs": list(
                    alignment.get_aligned_pairs(matches_only=True)
                ),
            }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--chromosome", required=True)
    parser.add_argument("--position", required=True, type=int)
    parser.add_argument("--ref", required=True)
    parser.add_argument("--alt", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--reverse-retain-bases",
        type=int,
        default=None,
        help="Exploratory reverse-read prefix length; not a validated QC rule.",
    )
    args = parser.parse_args()

    if args.output.exists():
        parser.error(f"Refusing to overwrite: {args.output}")

    if args.position < 1:
        parser.error("Position must be positive")

    ref = args.ref.upper()
    alt = args.alt.upper()

    if ref not in "ACGT" or alt not in "ACGT" or ref == alt:
        parser.error("REF and ALT must be distinct A/C/G/T bases")

    root = args.results_dir.expanduser().resolve()
    reference = args.reference.expanduser().resolve()

    with pysam.FastaFile(str(reference)) as fasta:
        if args.chromosome not in fasta.references:
            parser.error(
                f"Chromosome {args.chromosome} is absent from the reference"
            )

        if args.position > fasta.get_reference_length(args.chromosome):
            parser.error("Position exceeds chromosome length")

        actual_ref = fasta.fetch(
            args.chromosome, args.position - 1, args.position
        ).upper()

    if actual_ref != ref:
        parser.error(
            f"REF mismatch at {args.chromosome}:{args.position}: "
            f"provided {ref}, reference FASTA contains {actual_ref}"
        )

    print(
        f"Reference validation: PASS "
        f"({args.chromosome}:{args.position} {actual_ref})"
    )

    inventory = read_csv(root / f"{args.sample}_Inventory.csv")
    qc_rows = read_csv(
        root / "alignment" / f"{args.sample}_Generalized_Alignment_QC.csv"
    )

    qc = {row["Read_ID"]: row for row in qc_rows}

    if not inventory:
        raise ValueError("Inventory is empty")

    output_rows = []

    for row in inventory:
        if row["Sample"] != args.sample:
            raise ValueError("Inventory sample mismatch")

        read_id = row["Read_ID"]

        if read_id not in qc:
            raise ValueError(f"Missing alignment QC for {read_id}")

        record = SeqIO.read(row["AB1_path"], "abi")
        original = str(record.seq).upper()
        retained = original
        diagnostic_trim = False

        if (
            row["Direction"] == "R"
            and args.reverse_retain_bases is not None
        ):
            if not 1 <= args.reverse_retain_bases <= len(original):
                raise ValueError("Invalid reverse retained length")

            retained = original[:args.reverse_retain_bases]
            diagnostic_trim = True

        alignment = align_sequence(retained, read_id, reference)

        if alignment["chromosome"] != args.chromosome:
            raise ValueError(
                f"{read_id}: aligned to {alignment['chromosome']}, "
                f"not {args.chromosome}"
            )

        matches = [
            query_index
            for query_index, reference_index in alignment["aligned_pairs"]
            if reference_index == args.position - 1
        ]

        if len(matches) != 1:
            raise ValueError(
                f"{read_id}: expected one aligned base at target; "
                f"found {len(matches)}"
            )

        query_index = matches[0]

        original_index = (
            len(retained) - 1 - query_index
            if alignment["reverse"]
            else query_index
        )

        raw = record.annotations["abif_raw"]
        order = raw["FWO_1"]

        if isinstance(order, bytes):
            order = order.decode()

        if sorted(order) != sorted("ACGT"):
            raise ValueError(f"Unexpected dye channel order: {order}")

        traces = {
            base: np.asarray(raw[f"DATA{i}"], dtype=float)
            for i, base in zip(range(9, 13), order)
        }

        scan = int(raw["PLOC2"][original_index])

        trace_ref = ref.translate(COMPLEMENT) if alignment["reverse"] else ref
        trace_alt = alt.translate(COMPLEMENT) if alignment["reverse"] else alt

        signals = {
            base: float(trace[scan])
            for base, trace in traces.items()
        }

        ranked = sorted(signals.values(), reverse=True)
        primary = ranked[0]
        secondary = ranked[1]

        ref_signal = signals[trace_ref]
        alt_signal = signals[trace_alt]

        output_rows.append({
            "Sample": args.sample,
            "Read_ID": read_id,
            "Gene": row["Gene"],
            "Target": row["Target"],
            "Direction": row["Direction"],
            "Chromosome": args.chromosome,
            "Position_GRCh38_1based": args.position,
            "REF": ref,
            "ALT": alt,
            "AB1_base_1based": original_index + 1,
            "AB1_called_base": original[original_index],
            "AB1_base_quality": record.letter_annotations[
                "phred_quality"
            ][original_index],
            "Reference_channel": trace_ref,
            "Alternative_channel": trace_alt,
            "Reference_signal": round(ref_signal, 3),
            "Alternative_signal": round(alt_signal, 3),
            "ALT_fraction_of_REF_ALT_signal": round(
                alt_signal / (ref_signal + alt_signal), 6
            ) if ref_signal + alt_signal else "",
            "Secondary_primary_ratio": round(
                secondary / primary, 6
            ) if primary else "",
            "Alignment_QC_flag": qc[read_id]["Alignment_QC_flag"],
            "Diagnostic_trim_applied": diagnostic_trim,
            "Diagnostic_retained_bases": len(retained),
            "Diagnostic_alignment_MAPQ": alignment["mapq"],
            "Diagnostic_alignment_CIGAR": alignment["cigar"],
            "Validation_status": "UNVALIDATED",
        })

    columns = list(output_rows[0])

    args.output.parent.mkdir(parents=True, exist_ok=True)

    with args.output.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(output_rows)

    print("\nTRACE-AWARE SNV EVIDENCE — EXPLORATORY")
    print("=" * 70)
    print("Sample:", args.sample)
    print("Target:", f"{args.chromosome}:{args.position} {ref}>{alt}")

    for row in output_rows:
        print(
            row["Read_ID"],
            row["Direction"],
            "REF=", row["Reference_signal"],
            "ALT=", row["Alternative_signal"],
            "ALT_fraction=", row["ALT_fraction_of_REF_ALT_signal"],
            "Original_QC=", row["Alignment_QC_flag"],
            "Diagnostic_trim=", row["Diagnostic_trim_applied"],
        )

    print("\nSaved:", args.output)
    print("All observations are UNVALIDATED.")
    print("No genotype or pathogenicity classification performed.")


if __name__ == "__main__":
    main()
