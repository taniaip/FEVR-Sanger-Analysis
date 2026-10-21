#!/usr/bin/env python3
"""Multi-gene Sanger alignment QC. Not a clinical variant caller."""

import argparse
import csv
import re
import shutil
import subprocess
from collections import Counter
from pathlib import Path

GENE_CHROMOSOMES = {
    "ADAM10": "chr15",
    "CDH5": "chr16",
    "FZD4": "chr11",
    "FZD7": "chr2",
    "IQGAP1": "chr15",
    "KIF11": "chr10",
    "LRP5": "chr11",
    "NDP": "chrX",
    "OLFM3": "chr1",
    "TSPAN12": "chr7",
    "ZNF408": "chr11",
}

FIELDS = [
    "Sample", "Read_ID", "Gene", "Target", "Direction",
    "Parse_status", "Expected_chromosome", "Chromosome",
    "Start_GRCh38_1based", "End_GRCh38_1based",
    "Alignment_strand", "MAPQ", "CIGAR",
    "Trimmed_length", "Aligned_query_bases", "Aligned_fraction",
    "Deletion_operations", "Insertion_operations",
    "Indel_size_category", "Alignment_QC_flag", "Review_flags",
]


def load_fastq(path):
    lengths = {}
    with path.open() as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().strip()
            plus = handle.readline()
            quality = handle.readline().strip()

            if not plus or not plus.startswith("+"):
                raise ValueError("Malformed FASTQ separator")
            if not header.startswith("@") or len(sequence) != len(quality):
                raise ValueError("Malformed FASTQ record")

            read_id = header[1:].split()[0]
            if read_id in lengths:
                raise ValueError(f"Duplicate FASTQ ID: {read_id}")
            lengths[read_id] = len(sequence)
    return lengths


def load_inventory(path, sample):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {
            "Sample", "Read_ID", "Gene", "Target",
            "Direction", "Sequencing_group", "Parse_status",
        }
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Inventory is missing required columns")

        rows = {}
        for row in reader:
            if row["Sample"] != sample:
                continue
            if row["Read_ID"] in rows:
                raise ValueError(f"Duplicate inventory ID: {row['Read_ID']}")
            rows[row["Read_ID"]] = row
    return rows


def parse_sam(path):
    records = {}
    with path.open() as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                raise ValueError("Malformed SAM record")

            read_id = fields[0]
            flag = int(fields[1])
            if flag & (0x100 | 0x800):
                continue
            if read_id in records:
                raise ValueError(f"Multiple primary alignments: {read_id}")

            if flag & 0x4:
                records[read_id] = None
                continue

            operations = re.findall(r"(\d+)([MIDNSHP=X])", fields[5])
            if not operations or "".join(n + op for n, op in operations) != fields[5]:
                raise ValueError(f"Invalid CIGAR for {read_id}")

            aligned = sum(int(n) for n, op in operations if op in "MI=X")
            reference_span = sum(int(n) for n, op in operations if op in "MDN=X")
            start = int(fields[3])

            records[read_id] = {
                "Chromosome": fields[2],
                "Start_GRCh38_1based": start,
                "End_GRCh38_1based": start + reference_span - 1,
                "Alignment_strand": "-" if flag & 16 else "+",
                "MAPQ": int(fields[4]),
                "CIGAR": fields[5],
                "Aligned_query_bases": aligned,
                "Deletion_operations": ";".join(
                    n for n, op in operations if op == "D"
                ),
                "Insertion_operations": ";".join(
                    n for n, op in operations if op == "I"
                ),
            }
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastq", required=True, type=Path)
    parser.add_argument("--inventory", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--reference-build", required=True, choices=["GRCh38"])
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--min-aligned-fraction", type=float, default=0.8)
    args = parser.parse_args()

    for program in ("minimap2", "samtools"):
        if not shutil.which(program):
            parser.error(f"{program} not found in PATH")

    if not 0 <= args.min_aligned_fraction <= 1:
        parser.error("Minimum aligned fraction must be between 0 and 1")
    if not 0 <= args.min_mapq <= 255:
        parser.error("Minimum MAPQ must be between 0 and 255")

    fastq = args.fastq.expanduser().resolve()
    inventory_path = args.inventory.expanduser().resolve()
    reference = args.reference.expanduser().resolve()
    output = args.output_dir.expanduser().resolve()

    for path in (fastq, inventory_path, reference):
        if not path.is_file():
            parser.error(f"Missing input file: {path}")

    lengths = load_fastq(fastq)
    inventory = load_inventory(inventory_path, args.sample)
    if not lengths:
        parser.error("FASTQ contains no reads")

    for read_id in lengths:
        if read_id not in inventory:
            parser.error(f"Read missing from inventory: {read_id}")
        if inventory[read_id]["Sequencing_group"] != "GENOMIC_FEVR":
            parser.error(f"Non-genomic read in FASTQ: {read_id}")

    sam = output / f"{args.sample}_generalized_genomic.sam"
    bam = output / f"{args.sample}_generalized_genomic.sorted.bam"
    bai = Path(str(bam) + ".bai")
    qc = output / f"{args.sample}_Generalized_Alignment_QC.csv"

    for path in (sam, bam, bai, qc):
        if path.exists():
            parser.error(f"Output already exists; refusing overwrite: {path}")

    output.mkdir(parents=True, exist_ok=True)

    print("Aligning genomic Sanger reads...", flush=True)
    with sam.open("w") as handle:
        subprocess.run(
            ["minimap2", "-a", "-x", "sr", "--secondary=no",
             str(reference), str(fastq)],
            stdout=handle,
            check=True,
        )

    subprocess.run(["samtools", "sort", "-o", str(bam), str(sam)], check=True)
    subprocess.run(["samtools", "index", str(bam)], check=True)
    alignments = parse_sam(sam)

    rows = []
    for read_id in sorted(lengths):
        metadata = inventory[read_id]
        alignment = alignments.get(read_id)
        if read_id not in alignments:
            raise RuntimeError(f"No SAM record for {read_id}")

        gene = metadata["Gene"].strip().upper()
        expected = GENE_CHROMOSOMES.get(gene, "")
        reviews = []

        if metadata["Parse_status"] != "OK":
            reviews.append("METADATA_REVIEW")
        if not expected:
            reviews.append("UNKNOWN_GENE_CHROMOSOME")

        if alignment is None:
            alignment = {}
            status = "UNMAPPED"
            fraction = 0
        else:
            fraction = alignment["Aligned_query_bases"] / lengths[read_id]
            if expected and alignment["Chromosome"] != expected:
                status = "WRONG_CHROMOSOME"
            elif alignment["MAPQ"] < args.min_mapq:
                status = "LOW_MAPQ"
            elif fraction < args.min_aligned_fraction:
                status = "LOW_ALIGNED_FRACTION"
            elif not expected:
                status = "UNVERIFIED_CHROMOSOME"
            else:
                status = "PASS"

        indel_lengths = [
            int(length)
            for field in ("Deletion_operations", "Insertion_operations")
            for length in alignment.get(field, "").split(";")
            if length
        ]

        if not indel_lengths:
            indel_category = "NONE"
        elif max(indel_lengths) >= 5:
            indel_category = "LARGE_5BP_PLUS"
            reviews.append("LARGE_INDEL_REVIEW")
        else:
            indel_category = "SMALL_1_4BP"

        rows.append({
            "Sample": args.sample,
            "Read_ID": read_id,
            "Gene": gene,
            "Target": metadata["Target"],
            "Direction": metadata["Direction"],
            "Parse_status": metadata["Parse_status"],
            "Expected_chromosome": expected,
            "Chromosome": alignment.get("Chromosome", ""),
            "Start_GRCh38_1based": alignment.get("Start_GRCh38_1based", ""),
            "End_GRCh38_1based": alignment.get("End_GRCh38_1based", ""),
            "Alignment_strand": alignment.get("Alignment_strand", ""),
            "MAPQ": alignment.get("MAPQ", 0),
            "CIGAR": alignment.get("CIGAR", ""),
            "Trimmed_length": lengths[read_id],
            "Aligned_query_bases": alignment.get("Aligned_query_bases", 0),
            "Aligned_fraction": round(fraction, 4),
            "Deletion_operations": alignment.get("Deletion_operations", ""),
            "Insertion_operations": alignment.get("Insertion_operations", ""),
            "Indel_size_category": indel_category,
            "Alignment_QC_flag": status,
            "Review_flags": ";".join(sorted(set(reviews))),
        })

    with qc.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    print("\nGENERALIZED ALIGNMENT QC")
    print("=" * 65)
    print("Sample:", args.sample)
    print("Reads:", len(rows))
    for status, count in sorted(Counter(
        row["Alignment_QC_flag"] for row in rows
    ).items()):
        print(f"  {status}: {count}")

    print("Reads with review flags:", sum(bool(r["Review_flags"]) for r in rows))
    print("\nSaved:", qc)
    print("Saved:", bam)


if __name__ == "__main__":
    main()
