#!/usr/bin/env python3
"""Step 914: report flag-based Q20 candidate quality classification."""
import argparse
from variant_analysis_utils import TYPES, load_csv, save_csv, quality_pass


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root', default='.', help='FEVR-Sanger-Analysis repository root')
    args = ap.parse_args()
    from pathlib import Path
    root = Path(args.root).expanduser().resolve()
    rows = load_csv(root / 'validation/FEVR_Consolidated_Candidates.csv')
    summary = []
    for typ in TYPES:
        subset = [r for r in rows if r['Variant_type'].strip().upper() == typ]
        passed = sum(quality_pass(r) for r in subset)
        summary.append(dict(Variant_type=typ, Initial_candidates=len(subset), Quality_pass=passed,
                            Quality_fail=len(subset)-passed, Quality_pass_percent=round(100*passed/len(subset), 2) if subset else 0,
                            Criterion='At least one Q20 ALT-supporting read (flag proxy)' if typ == 'SNV' else 'Both immediate flanks Q20 or higher (flag proxy)'))
    out = root / 'validation/FEVR_Variant_Quality_Filter_Summary.csv'
    save_csv(out, summary)
    print('STEP 914 — QUALITY FILTER SUMMARY')
    for r in summary:
        print(f'{r["Variant_type"]:<10} initial={r["Initial_candidates"]:>5} pass={r["Quality_pass"]:>5} fail={r["Quality_fail"]:>5} pass%={r["Quality_pass_percent"]:>6.2f}')
    print('TOTAL:', sum(r['Initial_candidates'] for r in summary), 'PASS:', sum(r['Quality_pass'] for r in summary))
    print('Saved:', out)

if __name__ == '__main__':
    main()
