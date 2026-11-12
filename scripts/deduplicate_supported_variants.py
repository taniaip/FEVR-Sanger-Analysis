#!/usr/bin/env python3
"""Step 917: aggregate identical Tier 1 variants across samples (Tier 1 only)."""
import argparse
from collections import defaultdict
from pathlib import Path
from variant_analysis_utils import load_csv, save_csv, ids, truth


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root', default='.')
    args = ap.parse_args()
    root = Path(args.root).expanduser().resolve()
    rows = load_csv(root / 'validation/FEVR_Ranked_Variant_Candidates.csv')
    tier1 = [r for r in rows if int(r['Evidence_tier']) == 1]
    groups = defaultdict(list)
    for r in tier1:
        groups[(r['Chromosome'], r['Position_GRCh38_1based'], r['REF'], r['ALT'], r['Variant_type'])].append(r)
    output = []
    for (chrom,pos,ref,alt,typ), records in groups.items():
        samples, genes, warnings = defaultdict(list), set(), set()
        for r in records:
            samples[r['Sample']].append(r)
            genes.update(x.strip() for x in r.get('Genes_and_targets','').split(';') if x.strip())
            if r.get('Alignment_review_flags','').strip():
                warnings.add(r['Alignment_review_flags'].strip())
        ordered = sorted(samples, key=lambda s: (0,int(s)) if s.isdigit() else (1,s))
        per_sample, bidir = {}, {}
        for sample in ordered:
            all_ids, f, rev = set(), set(), set()
            for r in samples[sample]:
                all_ids.update(ids(r.get('Supporting_read_IDs','')))
                f.update(ids(r.get('Forward_read_IDs','')))
                rev.update(ids(r.get('Reverse_read_IDs','')))
            per_sample[sample] = len(all_ids)
            bidir[sample] = bool(f and rev)
        output.append(dict(Genes_and_targets=';'.join(sorted(genes)), Chromosome=chrom,
                           Position_GRCh38_1based=pos, REF=ref, ALT=alt,
                           GRCh38_variant=f'{chrom}:{pos} {ref}>{alt}', Variant_type=typ,
                           Number_of_samples=len(ordered), Sample_IDs=';'.join(ordered),
                           Total_supporting_reads=sum(per_sample.values()),
                           Supporting_reads_by_sample=';'.join(f'{s}:{per_sample[s]}' for s in ordered),
                           Bidirectional_sample_count=sum(bidir.values()), Best_evidence_tier=1,
                           Alignment_review_flags=';'.join(sorted(warnings)), Validation_status='UNVALIDATED'))
    output.sort(key=lambda r:(-r['Number_of_samples'],-r['Total_supporting_reads'],r['Chromosome'],int(r['Position_GRCh38_1based'])))
    dest = root / 'validation'
    save_csv(dest / 'FEVR_Unique_Bidirectional_Variants.csv', output)
    save_csv(dest / 'FEVR_Recurrent_Bidirectional_Variants.csv', [r for r in output if r['Number_of_samples']>=2], fields=list(output[0]) if output else None)
    print('VARIANT DEDUPLICATION — DEDUPLICATED TIER 1')
    print('Tier 1 sample records:',len(tier1),'Unique variants:',len(output),'Recurrent variants:',sum(r['Number_of_samples']>=2 for r in output))
    for r in output:
        print(f'{r["GRCh38_variant"]:<34} samples={r["Number_of_samples"]} reads={r["Total_supporting_reads"]} {r["Supporting_reads_by_sample"]}')

if __name__ == '__main__':
    main()
