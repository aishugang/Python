#!/usr/bin/env python3
from __future__ import annotations
import csv, json
from pathlib import Path

TARGET = [
    {'event':1,'onset_h':151.64,'average_slip_m':0.0121,'seismic_moment_Nm':1.185e10,'moment_magnitude':0.6812,'static_stress_drop_Pa':2.0e6,'max_slip_velocity_mps':0.6228,'rupture_duration_s':0.4753},
    {'event':2,'onset_h':375.43,'average_slip_m':0.0106,'seismic_moment_Nm':1.035e10,'moment_magnitude':0.6433,'static_stress_drop_Pa':2.057e6,'max_slip_velocity_mps':0.8348,'rupture_duration_s':0.4446},
]

def relerr(x: float, y: float) -> float:
    return abs(x-y)/max(abs(y),1e-30)

def main() -> None:
    case = Path('work/Table5_case')
    path = case/'Table5_reproduced.csv'
    if not path.exists():
        raise SystemExit(f'Missing {path}')
    with path.open(newline='') as f:
        rows=[{k:(int(v) if k=='event' else float(v)) for k,v in row.items()} for row in csv.DictReader(f)]
    report={'target':TARGET,'reproduced':rows,'comparison':[]}
    for i,t in enumerate(TARGET):
        if i>=len(rows):
            report['comparison'].append({'event':i+1,'status':'missing'})
            continue
        r=rows[i]
        metrics={}
        for key in t:
            if key!='event':
                metrics[key]={'target':t[key],'value':r[key],'relative_error':relerr(r[key],t[key])}
        report['comparison'].append({'event':i+1,'status':'present','metrics':metrics})
    (case/'Table5_comparison.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    lines=['# Table 5 reproduction comparison','',
           'The printed rupture-duration unit is treated as seconds: the reported sub-millisecond values are smaller than the paper minimum time step and conflict with Fig. 19.','']
    for item in report['comparison']:
        lines += [f"## Event {item['event']}"]
        if item['status']!='present':
            lines += ['Missing event.','']; continue
        lines += ['| Metric | Target | Reproduced | Relative error |','|---|---:|---:|---:|']
        for key,d in item['metrics'].items():
            lines.append(f"| {key} | {d['target']:.12g} | {d['value']:.12g} | {d['relative_error']:.6%} |")
        lines.append('')
    text='\n'.join(lines)
    (case/'Table5_comparison.md').write_text(text,encoding='utf-8')
    print(text)

if __name__=='__main__':
    main()
