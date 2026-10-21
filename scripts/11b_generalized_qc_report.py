#!/usr/bin/env python3
"""Generate descriptive QC reports for preliminary FEVR Sanger analysis."""

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def read_csv(path, sample):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"Missing CSV header: {path}")
        rows = list(reader)

    if any(row.get("Sample") != sample for row in rows):
        raise ValueError(f"Unexpected sample ID in {path}")

    return rows


def write_csv(path, columns, rows):
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite: {path}")

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    sample = args.sample
    base = args.results_dir.expanduser().resolve()
    output = args.output_dir.expanduser().resolve()

    paths = {
        "inventory": base / f"{sample}_Inventory.csv",
        "initial": base / "initial" / f"{sample}_Initial_Read_QC.csv",
        "trim": base / "trimmed" / f"{sample}_Quality_Trim_Summary.csv",
        "alignment": base / "alignment" / f"{sample}_Generalized_Alignment_QC.csv",
        "targets": base / "coverage" / f"{sample}_Generalized_Target_Summary.csv",
        "snvs": base / "variants" / f"{sample}_Preliminary_SNV_Read_Evidence.csv",
        "candidates": base / "candidates" / f"{sample}_Preliminary_Candidates.csv",
    }

    for label, path in paths.items():
        if not path.is_file():
            parser.error(f"Missing {label}: {path}")

    data = {
        label: read_csv(path, sample)
        for label, path in paths.items()
    }

    inventory = data["inventory"]
    initial = data["initial"]
    trimmed = data["trim"]
    alignment = data["alignment"]
    targets = data["targets"]
    snvs = data["snvs"]
    candidates = data["candidates"]

    def unique_ids(rows, label):
        ids = [row["Read_ID"] for row in rows]
        if len(ids) != len(set(ids)):
            raise ValueError(f"Duplicate Read_ID in {label}")
        return set(ids)

    inventory_ids = unique_ids(inventory, "inventory")
    initial_ids = unique_ids(initial, "initial QC")
    trim_ids = unique_ids(trimmed, "trimming QC")
    alignment_ids = unique_ids(alignment, "alignment QC")

    if inventory_ids != initial_ids or inventory_ids != trim_ids:
        raise ValueError("Inventory, initial QC and trimming IDs differ")

    if not alignment_ids.issubset(inventory_ids):
        raise ValueError("Alignment contains unknown read IDs")

    inventory_by_id = {r["Read_ID"]: r for r in inventory}

    for row in alignment:
        source = inventory_by_id[row["Read_ID"]]
        for field in ("Gene", "Target", "Parse_status"):
            if row[field] != source[field]:
                raise ValueError(
                    f"Alignment metadata mismatch: {row['Read_ID']} {field}"
                )

    passing = [
        r for r in alignment
        if r["Alignment_QC_flag"] == "PASS"
    ]
    passing_ids = {r["Read_ID"] for r in passing}

    if not passing_ids:
        raise ValueError(
            "INSUFFICIENT SEQUENCING EVIDENCE: no passing alignments"
        )

    for row in candidates:
        if row["Read_ID"] not in passing_ids:
            raise ValueError("Candidate from nonpassing read")
        if row["Validation_status"] != "UNVALIDATED":
            raise ValueError("Unexpected candidate validation status")

    for row in snvs:
        if row["Read_ID"] not in passing_ids:
            raise ValueError("SNV from nonpassing read")

    genes = sorted({
        r["Gene"] for r in inventory
        if r["Sequencing_group"] == "GENOMIC_FEVR"
    })

    gene_rows = []
    for gene in genes:
        gene_inventory = [
            r for r in inventory
            if r["Gene"] == gene
            and r["Sequencing_group"] == "GENOMIC_FEVR"
        ]
        gene_alignment = [r for r in alignment if r["Gene"] == gene]
        gene_passing = [r for r in passing if r["Gene"] == gene]
        gene_candidates = [r for r in candidates if r["Gene"] == gene]

        gene_rows.append({
            "Sample": sample,
            "Gene": gene,
            "Inventory_genomic_reads": len(gene_inventory),
            "Trimmed_genomic_reads": len(gene_alignment),
            "Passing_alignment_reads": len(gene_passing),
            "Passing_metadata_review_reads": sum(
                r["Parse_status"] != "OK" for r in gene_passing
            ),
            "Candidate_observations": len(gene_candidates),
            "Candidate_observations_from_metadata_review_reads": sum(
                r["Parse_status"] != "OK" for r in gene_candidates
            ),
            "Candidate_observations_with_alignment_review_flags": sum(
                bool(r["Alignment_review_flags"].strip())
                for r in gene_candidates
            ),
            "Evidence_status": (
                "PASSING_READS_PRESENT"
                if gene_passing else "NO_PASSING_READS"
            ),
            "Exon_completeness_assessed": "NO",
            "Variant_validation_status": "UNVALIDATED",
        })

    target_rows = []
    for row in targets:
        target_rows.append({
            "Sample": sample,
            "Gene": row["Gene"],
            "Target": row["Target"],
            "Total_trimmed_reads": row["Total_trimmed_reads"],
            "Passing_reads": row["Passing_reads"],
            "Covered_reference_positions": row["Covered_reference_positions"],
            "Maximum_read_depth": row["Maximum_read_depth"],
            "Mean_depth_at_covered_positions":
                row["Mean_depth_at_covered_positions"],
            "Evidence_status": (
                "PASSING_READS_PRESENT"
                if int(row["Passing_reads"]) > 0
                else "NO_PASSING_READS"
            ),
            "Exon_completeness_assessed":
                row["Exon_completeness_assessed"],
        })

    counts = Counter(r["Variant_type"] for r in candidates)

    sample_summary = [
        {"Metric": "Sample", "Value": sample},
        {"Metric": "Inventory_AB1_reads", "Value": len(inventory)},
        {"Metric": "Genomic_inventory_reads", "Value": sum(
            r["Sequencing_group"] == "GENOMIC_FEVR"
            for r in inventory
        )},
        {"Metric": "Metadata_review_inventory_reads", "Value": sum(
            r["Parse_status"] != "OK" for r in inventory
        )},
        {"Metric": "Trimmed_genomic_reads", "Value": len(alignment)},
        {"Metric": "Passing_alignment_reads", "Value": len(passing)},
        {"Metric": "Passing_metadata_review_reads", "Value": sum(
            r["Parse_status"] != "OK" for r in passing
        )},
        {"Metric": "Genes_in_genomic_inventory", "Value": len(genes)},
        {"Metric": "Genes_with_passing_reads", "Value": sum(
            r["Passing_alignment_reads"] > 0 for r in gene_rows
        )},
        {"Metric": "Genes_without_passing_reads", "Value": sum(
            r["Passing_alignment_reads"] == 0 for r in gene_rows
        )},
        {"Metric": "SNV_observations", "Value": counts["SNV"]},
        {"Metric": "Deletion_observations", "Value": counts["DELETION"]},
        {"Metric": "Insertion_observations", "Value": counts["INSERTION"]},
        {"Metric": "Total_candidate_observations", "Value": len(candidates)},
        {"Metric": "SNV_observations_below_Q20", "Value": sum(
            int(r["Base_quality"]) < 20 for r in snvs
        )},
        {"Metric": "Candidates_from_metadata_review_reads", "Value": sum(
            r["Parse_status"] != "OK" for r in candidates
        )},
        {"Metric": "Candidates_with_alignment_review_flags", "Value": sum(
            bool(r["Alignment_review_flags"].strip())
            for r in candidates
        )},
        {"Metric": "Exon_completeness_assessed", "Value": "NO"},
        {"Metric": "Variant_validation_status", "Value": "UNVALIDATED"},
        {"Metric": "Genotype_or_clinical_interpretation", "Value": "NOT_PERFORMED"},
    ]

    files = [
        (
            output / f"{sample}_Sample_QC_Summary.csv",
            ["Metric", "Value"],
            sample_summary,
        ),
        (
            output / f"{sample}_Gene_QC_Summary.csv",
            list(gene_rows[0]) if gene_rows else [
                "Sample", "Gene", "Inventory_genomic_reads",
                "Trimmed_genomic_reads", "Passing_alignment_reads",
                "Passing_metadata_review_reads",
                "Candidate_observations",
                "Candidate_observations_from_metadata_review_reads",
                "Candidate_observations_with_alignment_review_flags",
                "Evidence_status", "Exon_completeness_assessed",
                "Variant_validation_status",
            ],
            gene_rows,
        ),
        (
            output / f"{sample}_Target_QC_Summary.csv",
            list(target_rows[0]) if target_rows else [
                "Sample", "Gene", "Target", "Total_trimmed_reads",
                "Passing_reads", "Covered_reference_positions",
                "Maximum_read_depth", "Mean_depth_at_covered_positions",
                "Evidence_status", "Exon_completeness_assessed",
            ],
            target_rows,
        ),
    ]

    for path, _, _ in files:
        if path.exists():
            raise FileExistsError(f"Refusing to overwrite: {path}")

    output.mkdir(parents=True, exist_ok=True)

    for path, columns, rows in files:
        write_csv(path, columns, rows)

    print("\nGENERALIZED FEVR SAMPLE QC REPORT")
    print("=" * 65)
    print("Sample:", sample)
    print("Passing alignments:", len(passing))
    print("Genes with passing reads:",
          sum(r["Passing_alignment_reads"] > 0 for r in gene_rows))
    print("Genes without passing reads:",
          sum(r["Passing_alignment_reads"] == 0 for r in gene_rows))
    print("Preliminary candidate observations:", len(candidates))
    print("All observations remain UNVALIDATED.")
    print("\nSaved:")
    for path, _, _ in files:
        print(" ", path)


if __name__ == "__main__":
    main()
