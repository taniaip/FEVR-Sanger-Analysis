"""Shared helpers for FEVR Sanger evidence reporting."""
import csv
import re
from pathlib import Path

TYPES = ('SNV', 'DELETION', 'INSERTION')
KEYS = ('Sample', 'Variant_type', 'Chromosome', 'Position_GRCh38_1based', 'REF', 'ALT')
READ_RE = re.compile(r'\bread_\d+\b')


def load_csv(path):
    with Path(path).open(newline='', encoding='utf-8-sig') as fh:
        return list(csv.DictReader(fh))


def save_csv(path, rows, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows and not fields:
        raise ValueError(f'Cannot write empty CSV without columns: {path}')
    with path.open('w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=fields or list(rows[0]), extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def ids(value):
    return set(READ_RE.findall(value or ''))


def truth(value):
    return str(value).strip().lower() in ('true', '1', 'yes')


def key(row):
    return tuple(str(row.get(field, '')).strip() for field in KEYS)


def quality_pass(row):
    """Reproduce Step 914 FLAG-BASED proxy, not direct Phred re-evaluation."""
    flags = row.get('Evidence_flags', '')
    if row['Variant_type'] == 'SNV':
        return 'NO_QUALIFYING_ALT_SUPPORT' not in flags
    if row['Variant_type'] in ('DELETION', 'INSERTION'):
        return 'LOW_FLANK_QUALITY' not in flags
    raise ValueError(f'Unexpected variant type: {row["Variant_type"]}')


def root_from_args(args):
    return Path(args.root).expanduser().resolve()
