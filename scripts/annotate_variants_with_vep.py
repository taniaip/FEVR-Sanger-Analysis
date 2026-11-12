#!/usr/bin/env python3
"""Step 918: Ensembl VEP GRCh38 batch annotation; raw response retained.

Colocated ClinVar/gnomAD records are NOT assumed to match the exact alternate allele.
This script intentionally omits unverified allele-specific clinical/frequency claims.
"""
import argparse
import csv
import json
import time
import urllib.error
import urllib.request
from pathlib import Path
from variant_analysis_utils import load_csv, save_csv

API = 'https://rest.ensembl.org/vep/human/region'


def variant_input(r):
    chrom = r['Chromosome'].removeprefix('chr')
    return f'{chrom} {int(r["Position_GRCh38_1based"])} . {r["REF"].upper()} {r["ALT"].upper()} . . .'


def request_batch(variants, retries=5):
    url = API + '?canonical=1&mane=1&hgvs=1&variant_class=1&numbers=1&protein=1&sift=b&polyphen=b&af=1&af_gnomade=1&clin_sig=1'
    request = urllib.request.Request(url, data=json.dumps({'variants': variants}).encode(),
        headers={'Content-Type':'application/json','Accept':'application/json','User-Agent':'FEVR-Sanger-Analysis/1.0'},method='POST')
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(request,timeout=90) as response:
                data = json.load(response)
                if not isinstance(data,list):
                    raise ValueError('Expected VEP JSON array')
                return data
        except urllib.error.HTTPError as e:
            detail = e.read().decode(errors='replace')[:500]
            if e.code not in (429,500,502,503,504) or attempt == retries-1:
                raise RuntimeError(f'VEP HTTP {e.code}: {detail}') from e
            wait = min(60,3*2**attempt)
            print(f'HTTP {e.code}; retrying after {wait}s'); time.sleep(wait)
        except (TimeoutError,urllib.error.URLError):
            if attempt == retries-1:
                raise
            time.sleep(min(60,3*2**attempt))
    raise RuntimeError('VEP request failed')


def select_transcript(result, genes_and_targets):
    transcripts = result.get('transcript_consequences',[])
    genes = {x.split(':',1)[0].strip().upper() for x in genes_and_targets.split(';') if x.strip()}
    matching = [t for t in transcripts if t.get('gene_symbol','').upper() in genes]
    choices = matching or transcripts
    if not choices:
        return {}
    return next((t for t in choices if t.get('mane_select')),
           next((t for t in choices if t.get('canonical') == 1), choices[0]))


def stringify(value):
    return ';'.join(map(str,value)) if isinstance(value,list) else '' if value is None else str(value)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',default='.')
    ap.add_argument('--batch-size',type=int,default=20)
    ap.add_argument('--offline-raw',help='Reuse an existing VEP raw JSON file; no network requests')
    args = ap.parse_args()
    if args.batch_size < 1 or args.batch_size > 200:
        ap.error('--batch-size must be between 1 and 200')
    root = Path(args.root).expanduser().resolve()
    dest = root / 'validation'
    rows = load_csv(dest / 'FEVR_Unique_Bidirectional_Variants.csv')
    if not rows:
        raise SystemExit('No Tier 1 variants to annotate')
    inputs = [variant_input(r) for r in rows]
    if len(set(inputs)) != len(inputs):
        raise SystemExit('Duplicate VEP input strings; inspect normalized variants')
    if args.offline_raw:
        results = json.loads(Path(args.offline_raw).read_text(encoding='utf-8'))
    else:
        results = []
        for start in range(0,len(inputs),args.batch_size):
            batch = inputs[start:start+args.batch_size]
            print(f'Annotating {start+1}–{start+len(batch)} of {len(inputs)}')
            results.extend(request_batch(batch))
            time.sleep(1)
        # Only overwrite raw file after all requests succeeded.
        (dest / 'FEVR_Variant_VEP_Raw.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
    by_input = {}
    for r in results:
        inp = r.get('input','')
        if inp in by_input:
            raise SystemExit(f'Duplicate VEP response input: {inp}')
        by_input[inp] = r
    output = []
    for row, inp in zip(rows,inputs):
        result = by_input.get(inp,{})
        transcript = select_transcript(result,row.get('Genes_and_targets',''))
        colocated = result.get('colocated_variants',[])
        output.append({**row, 'VEP_input':inp, 'Most_severe_consequence':result.get('most_severe_consequence',''),
                       'MANE_transcript':stringify(transcript.get('mane_select')),
                       'Selected_transcript':transcript.get('transcript_id',''),
                       'HGVSc':transcript.get('hgvsc',''), 'HGVSp':transcript.get('hgvsp',''),
                       'Consequence_terms':stringify(transcript.get('consequence_terms',[])),
                       'SIFT':transcript.get('sift_prediction',''), 'PolyPhen':transcript.get('polyphen_prediction',''),
                       'Colocated_variant_IDs_UNVERIFIED':';'.join(sorted({str(v['id']) for v in colocated if v.get('id')})),
                       'gnomAD_frequency':'', 'ClinVar_significance':'',
                       'Population_clinical_annotation_status':'NOT_ALLELE_VERIFIED',
                       'VEP_annotation_found':bool(result)})
    out = dest / 'FEVR_Variant_VEP_Annotations.csv'
    save_csv(out,output)
    print('VEP ANNOTATION — VEP GRCh38')
    print('Input variants:',len(rows),'Annotated:',sum(r['VEP_annotation_found'] for r in output),'Missing:',sum(not r['VEP_annotation_found'] for r in output))
    for r in output:
        print(f'{r["GRCh38_variant"]:<35} {r["Most_severe_consequence"]:<36} {r["Selected_transcript"]}')
    print('Saved:',out)

if __name__ == '__main__':
    main()
