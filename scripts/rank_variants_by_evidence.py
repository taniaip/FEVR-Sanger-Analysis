#!/usr/bin/env python3
"""Step 916: rank Q20-passing sample-specific candidate records by read evidence."""
import argparse
from collections import Counter
from pathlib import Path
from variant_analysis_utils import TYPES, load_csv, save_csv, key, truth


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root', default='.')
    args = ap.parse_args()
    root = Path(args.root).expanduser().resolve()
    original_rows = load_csv(root / 'validation/FEVR_Consolidated_Candidates.csv')
    original = {key(r): r for r in original_rows}
    if len(original) != len(original_rows):
        raise SystemExit('Duplicate candidate keys in consolidated CSV: cannot safely merge by key.')
    audit = load_csv(root / 'validation/FEVR_Variant_Read_Support_Audit.csv')
    if len(audit) != len(original_rows):
        raise SystemExit('Audit row count differs from original candidate count.')
    hits = []
    for r in audit:
        if not truth(r['Q20_pass']):
            continue
        count = int(r['Supporting_read_count'])
        bidirectional = truth(r['Bidirectional_support'])
        tier = 1 if bidirectional else 2 if count >= 2 else 3 if count == 1 else 4
        base = original.get(key(r))
        if base is None:
            raise SystemExit(f'Unmatched audit key: {key(r)}')
        hits.append({**base, 'Evidence_tier': tier,
                     'Evidence_tier_description': {1:'Bidirectional F+R support', 2:'Multiple reads, no F+R support', 3:'Single supporting read', 4:'No supporting reads'}[tier],
                     'Q20_pass': r['Q20_pass'], 'Supporting_read_count': count,
                     'Forward_read_IDs': r['Forward_read_IDs'], 'Reverse_read_IDs': r['Reverse_read_IDs'],
                     'Unknown_direction_read_IDs': r['Unknown_direction_read_IDs'],
                     'Bidirectional_support': bidirectional, 'Evidence_IDs_match': r['Evidence_IDs_match']})
    hits.sort(key=lambda r: (r['Evidence_tier'], -r['Supporting_read_count'], r['Sample'],
                             r.get('Genes_and_targets',''), r['Chromosome'], int(r['Position_GRCh38_1based'])))
    summary = []
    counts = Counter((r['Evidence_tier'], r['Variant_type']) for r in hits)
    for tier in (1,2,3,4):
        entry = {'Tier':tier, **{typ:counts[(tier,typ)] for typ in TYPES}}
        entry['Total'] = sum(entry[typ] for typ in TYPES)
        summary.append(entry)
    dest = root / 'validation'
    save_csv(dest / 'FEVR_Ranked_Variant_Candidates.csv', hits)
    save_csv(dest / 'FEVR_Bidirectionally_Supported_Variants.csv', [r for r in hits if r['Evidence_tier']==1], fields=list(hits[0]))
    save_csv(dest / 'FEVR_Variant_Evidence_Tier_Summary.csv', summary)
    print('EVIDENCE RANKING — EVIDENCE TIERS')
    for r in summary:
        print(r)
    print('TOTAL RETAINED:', len(hits))

if __name__ == '__main__':
    main()
