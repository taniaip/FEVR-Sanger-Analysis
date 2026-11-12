#!/usr/bin/env python3
"""Run Steps 914–918 sequentially. By default stops before network-dependent VEP."""
import argparse
import subprocess
import sys
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',default='.')
    ap.add_argument('--with-vep',action='store_true',help='Run VEP REST annotation (requires internet)')
    ap.add_argument('--strict',action='store_true',help='Enforce strict read/evidence audit in Step 915')
    args = ap.parse_args()
    root = str(Path(args.root).expanduser().resolve())
    scripts = Path(__file__).resolve().parent
    names = ['filter_variants_by_quality.py','audit_variant_read_support.py','rank_variants_by_evidence.py','deduplicate_supported_variants.py']
    if args.with_vep:
        names.append('annotate_variants_with_vep.py')
    for name in names:
        command = [sys.executable,str(scripts/name),'--root',root]
        if name == 'audit_variant_read_support.py' and args.strict:
            command.append('--strict')
        print('\nRUNNING:',name,flush=True)
        subprocess.run(command,check=True)
    if not args.with_vep:
        print('\nSteps 914–917 complete. Re-run with --with-vep for online Step 918.')

if __name__ == '__main__':
    main()
