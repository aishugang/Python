#!/usr/bin/env python3
"""Score a completed calibration variant against the two published Table 5 events."""
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

# Timing and event count dominate this first-stage inverse calibration.
WEIGHTS = dict(onset_h=4.0, average_slip_m=1.0, seismic_moment_Nm=0.5,
               moment_magnitude=0.25, static_stress_drop_Pa=1.0,
               max_slip_velocity_mps=0.5, rupture_duration_s=0.25)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-dir", type=Path, required=True)
    args = parser.parse_args()
    csv_path = args.case_dir / "Table5_reproduced.csv"
    rows = []
    if csv_path.exists():
        with csv_path.open(newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                rows.append({k: float(v) for k, v in row.items() if k != "event"})

    score = 0.0
    details = []
    for i, target in enumerate(TARGETS):
        if i >= len(rows):
            score += 100.0
            details.append({"event": i + 1, "missing": True, "penalty": 100.0})
            continue
        actual = rows[i]
        terms = {}
        for key, target_value in target.items():
            rel = abs(actual[key] - target_value) / max(abs(target_value), 1e-30)
            weighted = WEIGHTS[key] * rel
            score += weighted
            terms[key] = {"target": target_value, "actual": actual[key],
                          "relative_error": rel, "weighted_error": weighted}
        details.append({"event": i + 1, "missing": False, "terms": terms})

    # Additional events are also inconsistent with the printed table.
    if len(rows) > 2:
        score += 10.0 * (len(rows) - 2)

    variant_path = args.case_dir / "VARIANT.json"
    variant = json.loads(variant_path.read_text(encoding="utf-8")) if variant_path.exists() else {}
    result = {"score": score, "event_count": len(rows), "variant": variant, "details": details}
    (args.case_dir / "CALIBRATION_SCORE.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"variant": variant.get("variant"), "event_count": len(rows), "score": score}, sort_keys=True))


if __name__ == "__main__":
    main()
