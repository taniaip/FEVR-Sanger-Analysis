# FEVR Sanger candidate evidence pipeline — Steps 914–918

Reproducible reporting scripts for the existing FEVR Sanger analysis. **No clinical validation is implied.** Input files are read without modification. The scripts write new or replace *generated* variant analysis workflow CSV/JSON outputs under `validation/`; keep previous versions under version control or back them up before rerunning.

## Prerequisites

- Python 3.9+ (standard library only)
- Repository containing `validation/FEVR_Consolidated_Candidates.csv`
- Per-sample inventory and evidence audits at `results/<sample>/Analysis/...`
- Internet connection for online Ensembl VEP (Step 918 only)

Copy `scripts/` into your repository (or keep it in a subdirectory). From the repository root:

```bash
python3 scripts/run_variant_analysis.py --root . --strict
python3 scripts/annotate_variants_with_vep.py --root .
```

Or run everything in one command:

```bash
python3 scripts/run_variant_analysis.py --root . --strict --with-vep
```

Run a single step with `python3 scripts/filter_variants_by_quality.py --root .` (similarly for steps 915–918). If you already have the VEP raw response, run Step 918 without network:

```bash
python3 scripts/annotate_variants_with_vep.py --root . --offline-raw validation/FEVR_Variant_VEP_Raw.json
```

## Output files

| Step | Output |
|---|---|
| 914 | `validation/FEVR_Variant_Quality_Filter_Summary.csv` |
| 915 | `validation/FEVR_Variant_Read_Support_Audit.csv`, `validation/FEVR_Variant_Read_Support_Summary.csv` |
| 916 | `validation/FEVR_Ranked_Variant_Candidates.csv`, `validation/FEVR_Bidirectionally_Supported_Variants.csv`, `validation/FEVR_Variant_Evidence_Tier_Summary.csv` |
| 917 | `validation/FEVR_Unique_Bidirectional_Variants.csv`, `validation/FEVR_Recurrent_Bidirectional_Variants.csv` |
| 918 | `validation/FEVR_Variant_VEP_Annotations.csv`, `validation/FEVR_Variant_VEP_Raw.json` (online mode) |

## Definitions and cautions

- **Step 914 reproduces the previously used *flag-based* Q20 classification**, not raw Phred recalculation: SNVs fail if `NO_QUALIFYING_ALT_SUPPORT` is present; indels fail if `LOW_FLANK_QUALITY` is present. SNVs use ALT-base quality, while indels use immediately adjacent flank quality. Check the original command-line parameters if the upstream evidence-generation thresholds were overridden.
- **Step 915** maps supporting read IDs to forward/reverse **primer directions** from the sample inventory. It cross-checks against exact-support IDs in the per-sample audits. Use `--strict` to halt on missing evidence files, missing direction mappings, or mismatches. Different samples may reuse `read_XXX` identifiers.
- **Step 916** assigns Tier 1 = F+R support, Tier 2 = >=2 supporting reads without F+R, Tier 3 = one supporting read, Tier 4 = none. These are *technical review tiers*, not ACMG/AMP classifications.
- **Step 917** aggregates exact `(chromosome, position, REF, ALT, type)` strings among **Tier 1 only**, counting distinct read IDs per sample. It does **not** normalize indels or include Tier 2/3 carriers. Genomic representation should be normalized before clinical cross-database comparisons.
- **Step 918** sends VCF-like alleles to Ensembl VEP `/vep/human/region` (GRCh38 endpoint), saves the raw response, and selects a gene-matched MANE transcript where available (otherwise canonical or first). `Most_severe_consequence` is the overall VEP consequence, potentially from a different transcript than `Selected_transcript`. Colocated IDs are explicitly **UNVERIFIED**; gnomAD and ClinVar fields remain blank pending **exact-allele verification**. The KIF11 exon 13 insertion's apparent frameshift consequence requires separate transcript/normalization investigation.
- `Validation_status` remains `UNVALIDATED`; bidirectional support does not prove independent PCR confirmation, true heterozygosity, pathogenicity, or segregation.

## Expected historical counts (from project outputs)

- 2,538 candidate records: 1,567 SNVs, 804 deletions, 167 insertions.
- 582 preliminary Q20 passes: 346 SNVs, 191 deletions, 45 insertions.
- 62 candidates with >=2 supporting reads; 43 with F+R support.
- Tier 1: 43 sample-specific records representing 23 exact variant strings.

Counts should be checked on every rerun; they are not hard-coded into calculations.

## GitHub hygiene

Do not commit original patient/sample identifiers, AB1/PHD chromatograms, personal filesystem paths, or potentially sensitive patient-level result CSVs to a public repository without the appropriate permissions and de-identification. Consider committing scripts and documentation only, with generated results ignored by `.gitignore`.
