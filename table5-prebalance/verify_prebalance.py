#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

TARGETS = [
    dict(onset_h=151.64, average_slip_m=0.0121, seismic_moment_Nm=1.185e10,
         moment_magnitude=0.6812, static_stress_drop_Pa=2.0e6,
         max_slip_velocity_mps=0.6228, rupture_duration_s=0.4753),
    dict(onset_h=375.43, average_slip_m=0.0106, seismic_moment_Nm=1.035e10,
         moment_magnitude=0.6433, static_stress_drop_Pa=2.057e6,
         max_slip_velocity_mps=0.8348, rupture_duration_s=0.4446),
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--case-dir', type=Path, required=True)
    args = parser.parse_args()
    path = args.case_dir / 'Table5_reproduced.csv'
    if not path.exists():
        raise SystemExit(f'Missing {path}')

    with path.open(newline='', encoding='utf-8') as f:
        rows = [{k: float(v) for k, v in row.items() if k != 'event'} for row in csv.DictReader(f)]

    comparison = []
    for i, target in enumerate(TARGETS):
        if i >= len(rows):
            comparison.append({'event': i + 1, 'status': 'missing'})
            continue
        actual = rows[i]
        metrics = {}
        for key, value in target.items():
            rel = abs(actual[key] - value) / max(abs(value), 1e-30)
            metrics[key] = {'target': value, 'actual': actual[key], 'relative_error': rel}
        comparison.append({'event': i + 1, 'status': 'present', 'metrics': metrics})

    report = {'event_count': len(rows), 'targets': TARGETS, 'reproduced': rows, 'comparison': comparison}
    (args.case_dir / 'Table5_comparison.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')

    lines = ['# Table 5 prebalance-reset comparison', '', f'Event count: {len(rows)}', '']
    for item in comparison:
        lines.append(f"## Event {item['event']}")
        if item['status'] == 'missing':
            lines += ['Missing event.', '']
            continue
        lines += ['| Metric | Target | Prebalanced result | Relative error |', '|---|---:|---:|---:|']
        for key, data in item['metrics'].items():
            lines.append(f"| {key} | {data['target']:.12g} | {data['actual']:.12g} | {data['relative_error']:.6%} |")
        lines.append('')
    text = '\n'.join(lines)
    (args.case_dir / 'Table5_comparison.md').write_text(text, encoding='utf-8')
    print(text)


if __name__ == '__main__':
    main()
