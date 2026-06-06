#!/usr/bin/env python3
"""Compute summary statistics over the final verification CSVs.
Reads results/J{1,2}_OmegaMid.csv and writes results/verification_summary.json and results/verification_summary.md.

All J_lower values are reported with full 17-digit precision (the format used
in the source CSVs), preserving the rigorous numerical content so that the
summary statistics are reproducible from the raw data.
"""
import csv, json, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def load(conj):
    with open(os.path.join(ROOT, 'results', conj+'_OmegaMid.csv')) as f:
        return list(csv.DictReader(f))

def safe_float(s):
    try: return float(s)
    except: return None

summary = {}
for conj in ['J1', 'J2']:
    rows = load(conj)
    total = len(rows)
    ver = sum(1 for r in rows if r['verified'] == '1')
    # Keep original string representations keyed by value for full-precision recovery.
    val_str = []
    for r in rows:
        v = safe_float(r['J_lower'])
        if v is None or v == float('-inf') or v == float('inf'):
            continue
        val_str.append((v, r['J_lower'].strip()))
    val_str.sort(key=lambda t: t[0])
    n = len(val_str)

    def pct(p):
        # p in [0, 1]; returns (float, original_string) at that percentile
        if n == 0: return (None, None)
        idx = min(n - 1, max(0, int(p * n)))
        return val_str[idx]

    pmin = val_str[0]
    pmax = val_str[-1]
    p1   = pct(0.01)
    p10  = pct(0.10)
    p50  = pct(0.50)
    p90  = pct(0.90)
    p99  = pct(0.99)

    stat = {
        'total_cells':    total,
        'verified_cells': ver,
        'verified_rate':  ver / total if total else 0,
        'J_lower': {
            'min':        pmin[1],
            'max':        pmax[1],
            'p1':         p1[1],
            'p10':        p10[1],
            'median':     p50[1],
            'p90':        p90[1],
            'p99':        p99[1],
        },
    }
    summary[conj] = stat

summary['overall'] = {
    'all_verified':        all(summary[c]['verified_cells'] == summary[c]['total_cells'] for c in ['J1', 'J2']),
    'cell_def_file':       'inputs/cell_def.csv',
}

with open(os.path.join(ROOT, 'results', 'verification_summary.json'), 'w') as f:
    json.dump(summary, f, indent=2)

lines = ['# Verification Summary', '']
lines.append('## Overall')
lines.append(f'- **Project fully verified:** {"YES" if summary["overall"]["all_verified"] else "NO"}')
lines.append(f'- Cell definition: `{summary["overall"]["cell_def_file"]}`')
lines.append('')
for conj in ['J1', 'J2']:
    s = summary[conj]
    ty = 'Laugesen–Siudeja type' if conj == 'J1' else 'Cheeger type'
    lines.append(f'## {conj} ({ty})')
    lines.append(f'- Total cells: **{s["total_cells"]:,}**')
    lines.append(f'- Verified:    **{s["verified_cells"]:,}**  ({s["verified_rate"] * 100:.3f}%)')
    lines.append('- `J_lower` statistics (full-precision values from the source CSV):')
    lines.append(f'  - **min**:    `{s["J_lower"]["min"]}`')
    lines.append(f'  - **median**: `{s["J_lower"]["median"]}`')
    lines.append(f'  - **max**:    `{s["J_lower"]["max"]}`')
    lines.append(f'  - 1%  percentile: `{s["J_lower"]["p1"]}`')
    lines.append(f'  - 10% percentile: `{s["J_lower"]["p10"]}`')
    lines.append(f'  - 90% percentile: `{s["J_lower"]["p90"]}`')
    lines.append(f'  - 99% percentile: `{s["J_lower"]["p99"]}`')
    lines.append('')
with open(os.path.join(ROOT, 'results', 'verification_summary.md'), 'w') as f:
    f.write('\n'.join(lines))
print(json.dumps(summary, indent=2))
