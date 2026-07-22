#!/usr/bin/env python3
"""Score fixed-Q runs using event count, first onset, and recurrence interval."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

T1_TARGET = 151.64
DT_TARGET = 375.43 - 151.64


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-dir", type=Path, required=True)
    args = parser.parse_args()

    rows = []
    csv_path = args.case_dir / "Table5_reproduced.csv"
    if csv_path.exists():
        with csv_path.open(newline="", encoding="utf-8") as f:
            rows = list(csv.DictReader(f))

    onset = [float(r["onset_h"]) for r in rows]
    score = 0.0
    terms = {}

    if not onset:
        score = 1000.0
        terms["missing_event1"] = 1000.0
    else:
        e1 = abs(onset[0] - T1_TARGET) / T1_TARGET
        score += 10.0 * e1
        terms["t1_relative_error"] = e1

    if len(onset) < 2:
        score += 250.0
        terms["missing_event2"] = 250.0
    else:
        recurrence = onset[1] - onset[0]
        edt = abs(recurrence - DT_TARGET) / DT_TARGET
        score += 10.0 * edt
        terms["recurrence_h"] = recurrence
        terms["recurrence_relative_error"] = edt

    if len(onset) > 2:
        extra = len(onset) - 2
        score += 50.0 * extra
        terms["extra_events_penalty"] = 50.0 * extra

    # Secondary Table 5 quantities are diagnostics only during this stage.
    secondary = []
    for row in rows[:2]:
        secondary.append({
            "average_slip_m": float(row["average_slip_m"]),
            "static_stress_drop_Pa": float(row["static_stress_drop_Pa"]),
            "max_slip_velocity_mps": float(row["max_slip_velocity_mps"]),
            "rupture_duration_s": float(row["rupture_duration_s"]),
        })

    variant_path = args.case_dir / "FIXEDQ_VARIANT.json"
    variant = json.loads(variant_path.read_text(encoding="utf-8")) if variant_path.exists() else {}
    result = {
        "score": score,
        "event_count": len(onset),
        "onset_h": onset,
        "target": {"t1_h": T1_TARGET, "recurrence_h": DT_TARGET},
        "terms": terms,
        "secondary": secondary,
        "variant": variant,
    }
    (args.case_dir / "FIXEDQ_SCORE.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"variant": variant.get("variant"), "event_count": len(onset), "score": score}, sort_keys=True))


if __name__ == "__main__":
    main()
