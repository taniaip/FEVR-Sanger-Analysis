# FEVR Sanger Sequencing Analysis

A research-use Python workflow for Sanger read quality control, genomic alignment, preliminary SNV/indel evidence assessment, and sample-, gene-, and target-level QC reporting in familial exudative vitreoretinopathy (FEVR) research.

> **Research use only.** This workflow does not provide validated genotypes, pathogenicity classifications, clinical diagnoses, or evidence of a variant-negative sample. All candidate observations remain `UNVALIDATED`.

## Generalized 11-stage workflow

Run `run_generalized_pipeline.py` from the repository root. The generalized pipeline uses the existing PHD extraction and quality trimming scripts, alongside the newer generalized analysis scripts.

| Stage | Script | Purpose |
| --- | --- | --- |
| 1 | `scripts/01b_generalized_inventory.py` | Inventory AB1/PHD pairs and parse gene, target, direction, and metadata review status. |
| 2 | `scripts/02_extract_phd_sequences.py` | Extract base calls and quality scores; report initial read QC. |
| 3 | `scripts/03_quality_trim.py` | Trim reads and separate genomic and plasmid sequencing groups. |
| 4 | `scripts/04b_generalized_alignment.py` | Align genomic reads with minimap2/samtools and assess alignment QC. |
| 5 | `scripts/05b_generalized_coverage.py` | Summarize passing reads and observed reference-position coverage. |
| 6 | `scripts/07b_generalized_variant_extraction.py` | Extract preliminary read-level SNV, deletion, and insertion observations. |
| 7 | `scripts/08b_generalized_snv_evidence.py` | Assess SNV quality and supporting/opposing read evidence. |
| 8 | `scripts/09b_generalized_deletion_evidence.py` | Assess deletion support and flanking base quality. |
| 9 | `scripts/09c_generalized_insertion_evidence.py` | Assess insertion support and flanking base quality. |
| 10 | `scripts/10b_generalized_candidates.py` | Consolidate read-level candidate observations and evidence flags. |
| 11 | `scripts/11b_generalized_qc_report.py` | Generate descriptive sample-, gene-, and target-level QC summaries. |

The original `run_pipeline.py` and numbered scripts are retained for the legacy workflow, including duplicate auditing, MAP3K6 plasmid primer QC, and Excel workbook generation.

## Requirements

- Python 3
- `minimap2` and `samtools` available on your `PATH`
- A suitable **independently verified** GRCh38 reference FASTA containing the chromosomes required for the targeted genes
- AB1 chromatograms (`.ab1`) and corresponding Phred files (`.phd.1`)

On macOS, users who do **not** already have the external tools can install them with:

```bash
brew install minimap2 samtools
```

If both tools are already installed, skip that step. Verify availability with:

```bash
minimap2 --version
samtools --version
```

The generalized runner uses the Python standard library plus these external command-line tools. The legacy Excel workbook workflow additionally uses `openpyxl`, listed in `requirements.txt`:

```bash
python3 -m pip install -r requirements.txt
```

## Run the generalized pipeline

From the repository root, with your own private sample and reference paths:

```bash
python3 run_generalized_pipeline.py \
  --sample SAMPLE_ID \
  --input /path/to/private/sample_directory \
  --reference /path/to/verified/GRCh38_reference.fa \
  --reference-build GRCh38 \
  --output-dir /path/to/private/results/SAMPLE_ID_generalized
```

The output directory must not exist before the run; the runner refuses to overwrite it. If filename metadata has been independently verified, an optional `--overrides /path/to/private/filename_overrides.csv` argument may be supplied. Do not guess metadata corrections.

The `--reference-build` label does **not** authenticate the reference sequence. Confirm its provenance, checksum, contig names, assembly, and target coordinates separately. The current analysis is designed around GRCh38 coordinates; changing the label is not a genome-build conversion.

## Gene and filename support

The generalized filename parser recognizes the labels `ADAM10`, `CDH5`, `FZD4`, `FZD7`, `GAPDH`, `IQGAP1`, `KIF11` (including `KIFF11`), `LRP5`, `MAP3K6`, `NDP`, `OLFM3`, `TSPAN12`, and `ZNF408` in the project's naming conventions. Recognition is not a guarantee of mapping, coverage, or biological identity.

Gene, target, and direction values derived from filenames are **parser assignments**, not verified laboratory metadata. Ambiguous records are marked `REVIEW`. Such reads may be retained in preliminary outputs; their uncertainty must be resolved before scientific interpretation.

The generalized runner analyzes genomic reads. The legacy MAP3K6 plasmid workflow is separate and does not establish plasmid variant calls. Some reused extraction/trimming outputs still have KIF11-specific filenames even though the files can contain genomic reads from other genes.

## Outputs and interpretation

The output directory contains `initial/`, `trimmed/`, `alignment/`, `coverage/`, `variants/`, `snv_evidence/`, `deletion_evidence/`, `insertion_evidence/`, `candidates/`, and `qc_report/`, plus the inventory CSV.

The `qc_report/` directory contains:

| File | Description |
| --- | --- |
| `<sample>_Sample_QC_Summary.csv` | Overall read counts, passing alignments, review flags, and candidate QC counts. |
| `<sample>_Gene_QC_Summary.csv` | Genomic read and candidate counts by gene, including genes without passing alignments. |
| `<sample>_Target_QC_Summary.csv` | Observed passing-read and reference-position coverage metrics by target. |

**`NO_PASSING_READS` means insufficient usable sequencing evidence, not a negative result.** The runner stops before downstream variant analysis if zero genomic reads pass alignment QC. If a sample has no observed SNVs, deletions, or insertions in a category, the corresponding evidence CSV is header-only rather than a pipeline failure.

Candidate CSV rows are **read-level observations**, not necessarily distinct variant alleles. Multiple read IDs do not prove biological independence. A passing alignment does not prove that each base call is reliable.

### Scientific limitations

- All candidate observations remain `UNVALIDATED`; no genotype or clinical interpretation is performed.
- Q20 thresholds and read-level evidence are preliminary screening criteria, not confirmation.
- Metadata marked `REVIEW` requires independent verification.
- Indels are not normalized; repeated-sequence contexts can produce ambiguous equivalent representations.
- Sequencing artifacts, low-quality bases, primer effects, and alignment ambiguity may generate apparent variants.
- Observed coverage does not establish complete exon or gene coverage; exon completeness is not assessed.
- No passing reads for a gene or target must never be interpreted as evidence of no variant.
- Independent confirmation and formal scientific/clinical validation remain outside this pipeline.

## Technical testing status

The generalized pipeline has completed internal research runs on three samples with different gene compositions and QC characteristics. Tests include exact reproducibility comparisons for the ten-stage pipeline and **17 CSV datasets** after integrating Stage 11, plus empty-category SNV/deletion/insertion handling and QC reporting for a gene with no passing reads.

These are **technical consistency tests**, not comprehensive biological validation. Synthetic fixtures, clean-environment testing, reference provenance checks, and additional independent scientific review are recommended before broader use.

## Data privacy and version control

Do not commit patient AB1/PHD files, sample metadata, identifying filenames, private paths, reference data without redistribution rights, or patient-level analysis results. The repository `.gitignore` excludes common sequencing and analysis formats and the `data/`, `results/`, and `references/` directories, but ignore rules are not a substitute for inspecting the staged files and their contents.

Before committing, review:

```bash
git diff --cached --name-only
git diff --cached
```

Use synthetic or appropriately authorized examples for any public repository.

## Legacy workflow

The original `run_pipeline.py` and scripts `01_inventory_ab1_phd.py` through `12_generate_workbook.py` remain available for the earlier research workflow. That workflow includes duplicate auditing, plasmid primer QC, and Excel workbook generation, and has distinct scope and output conventions from the generalized runner.

## License and release

No license has been selected. Confirm research-team ownership, institutional requirements, and release approval before making the repository public.
