# FEVR Sanger Sequencing Analysis

Python scripts for Sanger sequencing quality control, genomic alignment, and **preliminary, unvalidated** variant-evidence assessment in familial exudative vitreoretinopathy (FEVR) research.

> **Research use only.** This workflow does not provide validated genotypes, clinical classifications, or diagnostic results. All candidate variants require independent review and confirmation.

## Workflow

Run the scripts in numerical order, providing each script with the input files produced by the preceding step.

| Script | Purpose |
| --- | --- |
| `01_inventory_ab1_phd.py` | Inventory AB1 chromatograms and matching Phred PHD files; parse filename metadata. |
| `02_extract_phd_sequences.py` | Extract sequences and base qualities from PHD files and produce initial QC outputs. |
| `03_quality_trim.py` | Quality-trim sequencing reads and separate genomic and plasmid reads. |
| `04_align_genomic_reads.py` | Align genomic reads to a user-supplied reference using minimap2 and samtools. |
| `05_alignment_coverage.py` | Summarize passing alignments and reference-position coverage. |
| `06_duplicate_audit.py` | Audit AB1 files for byte-identical duplicates using SHA-256. |
| `07_extract_variants.py` | Extract preliminary SNV and indel observations from aligned reads. |
| `08_snv_evidence.py` | Assess base-quality and read-level reference/alternate evidence for SNVs. |
| `09_indel_evidence.py` | Assess preliminary deletion evidence and flanking base quality. |
| `10_consolidate_candidates.py` | Combine candidate observations and evidence summaries. |
| `11_plasmid_primer_qc.py` | Summarize MAP3K6 plasmid read quality by primer; does not call plasmid variants. |
| `12_generate_workbook.py` | Consolidate analysis tables into a multi-sheet Excel workbook. |

## Requirements

- Python 3
- Python package: `openpyxl>=3.1,<4` (listed in `requirements.txt`)
- External command-line tools: **minimap2** and **samtools**

Install the Python dependency from the repository root:

```bash
python3 -m pip install -r requirements.txt
```

Install minimap2 and samtools separately. For example, on macOS with Homebrew:

```bash
brew install minimap2 samtools
```

Check their availability with `minimap2 --version` and `samtools --version`.

## Inputs and filename conventions

The inventory stage expects Sanger chromatograms (`.ab1`) and corresponding Phred files (`.phd.1`). A suitable reference FASTA is required for genomic alignment.

The current inventory parser recognizes genomic labels `CDH5`, `KIF11`, `LRP5`, `NDP`, and `TSPAN12`, and MAP3K6 plasmid primer labels `F1`–`F6`, `R1`–`R6`, `M13F`, and `M13R`. Parsing depends on the filenames following the expected naming conventions; other naming schemes may require code changes.

### Example: inventory

```bash
python3 scripts/01_inventory_ab1_phd.py \
  --input /path/to/sample_directory \
  --sample SAMPLE_ID \
  --output /path/to/private_results/inventory.csv
```

The inventory CSV contains absolute paths and original filenames. Treat it as potentially sensitive.

For the remaining scripts, inspect their actual command-line options before running them:

```bash
python3 scripts/02_extract_phd_sequences.py --help
```

Repeat with the relevant script name. This README intentionally does not prescribe unverified command-line options for later stages. The scripts have been run sequentially for an initial research sample; a clean-environment end-to-end test is still needed before claiming general portability.

## Reference and analysis scope

The initial development run used a GRCh38 chromosome 10 reference for **KIF11** genomic reads. Confirm the reference assembly, contig identifiers, target coordinates, and strand conventions before interpreting any result. Some output labels are specific to KIF11 exon 10 and GRCh38; the workflow has **not** been demonstrated to be fully gene-agnostic.

The MAP3K6 plasmid branch performs primer-level quality control only; it does not compare plasmid sequences against a verified plasmid reference or call plasmid variants.

## Output interpretation

The workbook generator creates 13 worksheets: `Summary`, `Candidates`, `Candidate_Summary`, `SNV_Evidence`, `Indel_Quality`, `Raw_SNVs`, `Raw_Indels`, `Target_Coverage`, `Alignment_Details`, `Alignment_QC`, `Duplicate_Audit`, `Plasmid_Primer_QC`, and `Trim_QC`.

`Candidates` contains **read-level observations**, not necessarily unique alleles. Repeated observations of the same allele can occupy multiple rows.

### Scientific limitations

- **All preliminary variant observations remain `UNVALIDATED`.** They are not clinical findings.
- Different AB1 file hashes establish only that files are not byte-identical; they **do not establish biological independence** of the sequencing evidence.
- Multiple read IDs or forward/reverse primer labels do not, by themselves, establish independent biological confirmation.
- Base-quality thresholds and read-level evidence are screening criteria, not variant validation.
- Preliminary indel representations are **not normalized**; equivalent events may have different representations.
- The indel-quality assessment focuses on deletions; insertion validation is not implemented.
- Alignment artifacts, low-quality sequence, and primer-related effects may produce apparent differences.
- No pathogenicity classification or clinical interpretation is performed.

Independent confirmation, scientific review, and any necessary clinical validation are outside the scope of these scripts.

## Data privacy and version control

The repository `.gitignore` excludes raw chromatograms, PHD files, sequence files, alignments, reference FASTA files, CSV tables, Excel workbooks, and `data/`, `results/`, and `references/` directories. These exclusions help prevent accidental commits but **do not replace a manual review** of staged files or institutional data-sharing approval.

Before any public push, inspect staged files using `git diff --cached --name-only` and `git diff --cached`. Do not commit patient-identifiable metadata, private local paths, raw data, or sample-specific results without explicit authorization.

## Reproducibility status

The scripts were developed and executed for research Sample 426, with outputs retained locally and excluded from public version control. Reproducing the initial run requires appropriate access to the underlying research data and reference. Synthetic test fixtures, automated tests, a clean-environment run, and final scientific code review remain recommended before formal release.

## Code availability and citation

GitHub repository URL: **To be added after publication.**

Zenodo DOI: **To be added after archival.**

## License

**Not yet selected.** Confirm ownership and release permissions with the research team before choosing a license or making the repository public.


## Reference-build and target configuration

The genomic analysis pipeline currently supports GRCh38 coordinates only. Script 04 requires --reference-build GRCh38 alongside --reference. This argument declares the build but does not verify the FASTA sequence; users must verify the reference source and checksum independently.

Script 12 supports --gene, --target, and --reference-label. Its defaults preserve the Sample 426 KIF11 E10 workflow. Changing the reference label does not convert genomic coordinates or validate the reference sequence.
