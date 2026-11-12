#!/usr/bin/env python3
"""Step 915: audit supporting read IDs and forward/reverse primer directions."""
import argparse
from collections import Counter, defaultdict
from pathlib import Path
from variant_analysis_utils import TYPES, KEYS, load_csv, save_csv, ids, key, quality_pass

SPECS = {
    'SNV': ('snv_evidence', 'SNV_Evidence_Audit.csv', 'Supporting_read_IDs'),
    'DELETION': ('deletion_evidence', 'Indel_Quality_Evidence_Audit.csv', 'Exact_deletion_support_reads'),
    'INSERTION': ('insertion_evidence', 'Insertion_Evidence_Audit.csv', 'Exact_insertion_support_reads'),
}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root', default='.')
    ap.add_argument('--strict', action='store_true', help='Fail if evidence files, read directions or supporting IDs are inconsistent')
    args = ap.parse_args()
    root = Path(args.root).expanduser().resolve()
    rows = load_csv(root / 'validation/FEVR_Consolidated_Candidates.csv')
    samples = sorted({r['Sample'] for r in rows})
    directions, conflicts, missing = {}, set(), []
    for sample in samples:
        path = root / f'results/{sample}/Analysis/{sample}_Inventory.csv'
        if not path.exists():
            missing.append(str(path)); continue
        for r in load_csv(path):
            rid, direction = r.get('Read_ID', '').strip(), r.get('Direction', '').strip().upper()
            if not rid or direction not in ('F', 'R'):
                continue
            k = (sample, rid)
            if k in directions and directions[k] != direction:
                conflicts.add(k)
            else:
                directions[k] = direction
    for k in conflicts:
        directions.pop(k, None)
    evidence = defaultdict(set)
    evidence_keys = set()
    files_found = Counter()
    for sample in samples:
        for typ, (folder, suffix, column) in SPECS.items():
            path = root / f'results/{sample}/Analysis/{folder}/{sample}_{suffix}'
            if not path.exists():
                missing.append(str(path)); continue
            files_found[typ] += 1
            for r in load_csv(path):
                r['Sample'], r['Variant_type'] = sample, typ
                k = key(r)
                evidence_keys.add(k)
                evidence[k].update(ids(r.get(column, '')))
    audit = []
    for r in rows:
        typ = r['Variant_type'].strip().upper()
        if typ not in TYPES:
            raise ValueError(f'Unexpected variant type: {typ}')
        passed = quality_pass(r)
        support = ids(r.get('Supporting_read_IDs', ''))
        k = key(r)
        expected = evidence.get(k, set())
        forward = sorted(x for x in support if directions.get((r['Sample'], x)) == 'F')
        reverse = sorted(x for x in support if directions.get((r['Sample'], x)) == 'R')
        unknown = sorted(support - set(forward) - set(reverse))
        audit.append({**{field: r.get(field, '') for field in KEYS},
                      'Evidence_flags': r.get('Evidence_flags', ''), 'Q20_pass': passed,
                      'Supporting_read_IDs': ';'.join(sorted(support)), 'Supporting_read_count': len(support),
                      'Forward_read_IDs': ';'.join(forward), 'Reverse_read_IDs': ';'.join(reverse),
                      'Unknown_direction_read_IDs': ';'.join(unknown), 'Bidirectional_support': bool(forward and reverse),
                      'Audit_supporting_read_IDs': ';'.join(sorted(expected)),
                      'Evidence_audit_record_found': k in evidence_keys,
                      'Evidence_IDs_match': k in evidence_keys and support == expected})
    summary = []
    for typ in (*TYPES, 'TOTAL'):
        subset = [r for r in audit if typ == 'TOTAL' or r['Variant_type'] == typ]
        passed = [r for r in subset if r['Q20_pass']]
        summary.append(dict(Variant_type=typ, Initial_candidates=len(subset), Q20_pass=len(passed),
                            At_least_1_supporting_read=sum(r['Supporting_read_count'] >= 1 for r in passed),
                            At_least_2_supporting_reads=sum(r['Supporting_read_count'] >= 2 for r in passed),
                            Bidirectional_support=sum(r['Bidirectional_support'] for r in passed),
                            At_least_2_and_bidirectional=sum(r['Supporting_read_count'] >= 2 and r['Bidirectional_support'] for r in passed),
                            Unknown_direction_records=sum(bool(r['Unknown_direction_read_IDs']) for r in passed),
                            Evidence_mismatch_records=sum(not r['Evidence_IDs_match'] for r in passed)))
    save_csv(root / 'validation/FEVR_Variant_Read_Support_Audit.csv', audit)
    save_csv(root / 'validation/FEVR_Variant_Read_Support_Summary.csv', summary)
    print('STEP 915 — READ SUPPORT')
    for r in summary:
        print(f'{r["Variant_type"]:<10} initial={r["Initial_candidates"]} Q20={r["Q20_pass"]} >=2={r["At_least_2_supporting_reads"]} F+R={r["Bidirectional_support"]}')
    print('Evidence files:', dict(files_found), 'Missing files:', len(missing), 'Conflicting direction mappings:', len(conflicts))
    print('Q20-pass evidence mismatches:', summary[-1]['Evidence_mismatch_records'], 'Unknown direction records:', summary[-1]['Unknown_direction_records'])
    if args.strict and (missing or conflicts or summary[-1]['Evidence_mismatch_records'] or summary[-1]['Unknown_direction_records']):
        raise SystemExit('Strict audit failed; inspect outputs and missing files before proceeding.')

if __name__ == '__main__':
    main()
