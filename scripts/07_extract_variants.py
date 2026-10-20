#!/usr/bin/env python3

"""
Extract preliminary SNV and indel observations from passing Sanger reads.

Inputs:
    Sorted genomic BAM
    Alignment QC CSV
    Reference FASTA

Outputs:
    <sample>_Preliminary_SNV_Read_Evidence.csv
    <sample>_Preliminary_Indel_Read_Evidence.csv

Coordinates are 1-based, reference-oriented.
Indels are represented using adjacent reference anchoring where possible.
Indels are NOT normalized, independently confirmed, or clinically classified.

Requires samtools.
"""

import argparse
import csv
import re
import shutil
import subprocess
from pathlib import Path


def run_command(command):
    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout


def parse_cigar(cigar):
    matches = re.findall(r"(\d+)([MIDNSHP=X])", cigar)

    if not matches:
        raise ValueError(f"Invalid CIGAR: {cigar}")

    reconstructed = "".join(
        f"{length}{operation}"
        for length, operation in matches
    )

    if reconstructed != cigar:
        raise ValueError(f"Malformed CIGAR: {cigar}")

    return [(int(length), op) for length, op in matches]


def load_qc(path, sample):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))

    selected = {}

    for row in rows:
        if row["Sample"] != sample:
            continue

        if row["Alignment_QC_flag"] != "PASS":
            continue

        read_id = row["Read_ID"]

        if read_id in selected:
            raise ValueError(f"Duplicate QC read ID: {read_id}")

        selected[read_id] = row

    if not selected:
        raise ValueError("No passing genomic reads")

    return selected


def fetch_reference(reference, chromosome, start, end):
    """Fetch 1-based inclusive reference interval."""
    if start < 1 or end < start:
        raise ValueError("Invalid reference interval")

    region = f"{chromosome}:{start}-{end}"

    sequence = run_command([
        "samtools", "faidx", str(reference), region
    ])

    return "".join(
        line.strip()
        for line in sequence.splitlines()
        if not line.startswith(">")
    ).upper()


def write_csv(path, columns, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Extract raw Sanger SNV and indel observations."
    )

    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--alignment-qc", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)

    args = parser.parse_args()

    if shutil.which("samtools") is None:
        parser.error("samtools is not available")

    bam = args.bam.expanduser().resolve()
    qc_path = args.alignment_qc.expanduser().resolve()
    reference = args.reference.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()

    for path in [bam, qc_path, reference]:
        if not path.is_file():
            parser.error(f"Missing input: {path}")

    qc = load_qc(qc_path, args.sample)

    # Index the reference if required.
    if not Path(str(reference) + ".fai").is_file():
        subprocess.run(
            ["samtools", "faidx", str(reference)],
            check=True,
        )

    sam_text = run_command([
        "samtools", "view",
        "-F", "2308",
        str(bam),
    ])

    snvs = []
    indels = []
    seen = set()

    for line in sam_text.splitlines():
        fields = line.split("\t")

        if len(fields) < 11:
            raise ValueError("Malformed SAM record")

        read_id = fields[0]

        if read_id not in qc:
            continue

        if read_id in seen:
            raise ValueError(
                f"Multiple primary alignments for {read_id}"
            )

        seen.add(read_id)

        flag = int(fields[1])
        chromosome = fields[2]
        start = int(fields[3])
        mapq = int(fields[4])
        cigar = fields[5]
        sequence = fields[9].upper()
        quality_string = fields[10]

        if sequence == "*" or quality_string == "*":
            raise ValueError(
                f"Missing sequence or qualities: {read_id}"
            )

        qualities = [ord(char) - 33 for char in quality_string]

        if len(sequence) != len(qualities):
            raise ValueError(
                f"Sequence/quality mismatch: {read_id}"
            )

        operations = parse_cigar(cigar)

        reference_length = sum(
            length
            for length, op in operations
            if op in {"M", "D", "N", "=", "X"}
        )

        if reference_length < 1:
            continue

        reference_sequence = fetch_reference(
            reference,
            chromosome,
            start,
            start + reference_length - 1,
        )

        if len(reference_sequence) != reference_length:
            raise ValueError(
                f"Reference length mismatch for {read_id}"
            )

        strand = "-" if flag & 16 else "+"
        metadata = qc[read_id]

        query_index = 0
        reference_index = 0

        for length, operation in operations:

            if operation in {"M", "=", "X"}:
                for offset in range(length):
                    qpos = query_index + offset
                    rpos = reference_index + offset

                    ref_base = reference_sequence[rpos]
                    alt_base = sequence[qpos]

                    if (
                        ref_base in "ACGT"
                        and alt_base in "ACGT"
                        and ref_base != alt_base
                    ):
                        snvs.append({
                            "Sample": args.sample,
                            "Read_ID": read_id,
                            "Gene": metadata["Gene"],
                            "Target": metadata["Target"],
                            "Direction": metadata["Direction"],
                            "Chromosome": chromosome,
                            "Position_GRCh38_1based": start + rpos,
                            "REF": ref_base,
                            "ALT": alt_base,
                            "Variant_type": "SNV",
                            "Base_quality": qualities[qpos],
                            "Alignment_strand": strand,
                            "MAPQ": mapq,
                            "CIGAR": cigar,
                            "Validation_status": "UNVALIDATED",
                        })

                query_index += length
                reference_index += length

            elif operation == "I":
                inserted = sequence[
                    query_index:query_index + length
                ]

                # Left-anchor insertion when a preceding reference
                # base exists in the aligned reference interval.
                if reference_index > 0:
                    anchor = reference_sequence[reference_index - 1]
                    position = start + reference_index - 1

                    ref_allele = anchor
                    alt_allele = anchor + inserted
                    representation = "LEFT_ANCHORED_UNNORMALIZED"

                else:
                    # Fetch preceding genomic base if insertion occurs
                    # at the beginning of the alignment.
                    position = start - 1

                    if position < 1:
                        raise ValueError(
                            "Cannot anchor insertion at reference start"
                        )

                    anchor = fetch_reference(
                        reference, chromosome, position, position
                    )

                    ref_allele = anchor
                    alt_allele = anchor + inserted
                    representation = "LEFT_ANCHORED_UNNORMALIZED"

                left_quality = (
                    qualities[query_index - 1]
                    if query_index > 0 else ""
                )
                right_quality = (
                    qualities[query_index + length]
                    if query_index + length < len(qualities)
                    else ""
                )

                indels.append({
                    "Sample": args.sample,
                    "Read_ID": read_id,
                    "Gene": metadata["Gene"],
                    "Target": metadata["Target"],
                    "Direction": metadata["Direction"],
                    "Chromosome": chromosome,
                    "Position_GRCh38_1based": position,
                    "REF": ref_allele,
                    "ALT": alt_allele,
                    "Variant_type": "INSERTION",
                    "Indel_length": length,
                    "Query_index_0based": query_index,
                    "Read_length": len(sequence),
                    "Left_flank_quality": left_quality,
                    "Right_flank_quality": right_quality,
                    "Alignment_strand": strand,
                    "MAPQ": mapq,
                    "CIGAR": cigar,
                    "Representation": representation,
                    "Validation_status": "UNVALIDATED",
                })

                query_index += length

            elif operation == "D":
                deleted = reference_sequence[
                    reference_index:reference_index + length
                ]

                if reference_index > 0:
                    anchor = reference_sequence[reference_index - 1]
                    position = start + reference_index - 1

                else:
                    position = start - 1

                    if position < 1:
                        raise ValueError(
                            "Cannot anchor deletion at reference start"
                        )

                    anchor = fetch_reference(
                        reference, chromosome, position, position
                    )

                indels.append({
                    "Sample": args.sample,
                    "Read_ID": read_id,
                    "Gene": metadata["Gene"],
                    "Target": metadata["Target"],
                    "Direction": metadata["Direction"],
                    "Chromosome": chromosome,
                    "Position_GRCh38_1based": position,
                    "REF": anchor + deleted,
                    "ALT": anchor,
                    "Variant_type": "DELETION",
                    "Indel_length": length,
                    "Query_index_0based": query_index,
                    "Read_length": len(sequence),
                    "Left_flank_quality": (
                        qualities[query_index - 1]
                        if query_index > 0 else ""
                    ),
                    "Right_flank_quality": (
                        qualities[query_index]
                        if query_index < len(qualities)
                        else ""
                    ),
                    "Alignment_strand": strand,
                    "MAPQ": mapq,
                    "CIGAR": cigar,
                    "Representation": "LEFT_ANCHORED_UNNORMALIZED",
                    "Validation_status": "UNVALIDATED",
                })

                reference_index += length

            elif operation == "N":
                reference_index += length

            elif operation == "S":
                query_index += length

            elif operation in {"H", "P"}:
                pass

        if query_index != len(sequence):
            raise ValueError(
                f"CIGAR query length mismatch for {read_id}"
            )

        if reference_index != reference_length:
            raise ValueError(
                f"CIGAR reference length mismatch for {read_id}"
            )

    missing = set(qc) - seen

    if missing:
        raise RuntimeError(
            f"Passing reads missing from BAM: {sorted(missing)}"
        )

    output_dir.mkdir(parents=True, exist_ok=True)

    snv_path = output_dir / (
        f"{args.sample}_Preliminary_SNV_Read_Evidence.csv"
    )
    indel_path = output_dir / (
        f"{args.sample}_Preliminary_Indel_Read_Evidence.csv"
    )

    snv_columns = [
        "Sample", "Read_ID", "Gene", "Target", "Direction",
        "Chromosome", "Position_GRCh38_1based", "REF", "ALT",
        "Variant_type", "Base_quality", "Alignment_strand",
        "MAPQ", "CIGAR", "Validation_status"
    ]

    indel_columns = [
        "Sample", "Read_ID", "Gene", "Target", "Direction",
        "Chromosome", "Position_GRCh38_1based", "REF", "ALT",
        "Variant_type", "Indel_length", "Query_index_0based",
        "Read_length", "Left_flank_quality", "Right_flank_quality",
        "Alignment_strand", "MAPQ", "CIGAR", "Representation",
        "Validation_status"
    ]

    write_csv(snv_path, snv_columns, snvs)
    write_csv(indel_path, indel_columns, indels)

    print("\nPRELIMINARY SANGER VARIANT OBSERVATIONS")
    print("=" * 65)
    print("Sample:", args.sample)
    print("Passing genomic reads examined:", len(seen))
    print("Raw SNV observations:", len(snvs))
    print(
        "Distinct SNV alleles:",
        len({
            (r["Chromosome"], r["Position_GRCh38_1based"],
             r["REF"], r["ALT"])
            for r in snvs
        }),
    )
    print("Raw indel observations:", len(indels))
    print(
        "  Deletions:",
        sum(r["Variant_type"] == "DELETION" for r in indels)
    )
    print(
        "  Insertions:",
        sum(r["Variant_type"] == "INSERTION" for r in indels)
    )

    print("\nObservations by read:")
    for read_id in sorted(seen):
        snv_count = sum(r["Read_ID"] == read_id for r in snvs)
        indel_count = sum(r["Read_ID"] == read_id for r in indels)
        print(
            f"  {read_id}: {snv_count} SNVs, {indel_count} indels"
        )

    print("\nSaved:")
    print(" ", snv_path)
    print(" ", indel_path)
    print("\nAll observations are UNVALIDATED.")
    print("No genotype or pathogenicity classification performed.")


if __name__ == "__main__":
    main()
