#!/usr/bin/env python3
"""Apply a controlled parameter variant to a prepared Table 5 case."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

VARIANTS = {
    # Paper Table 4 hydraulic aperture, injection-rate timing sweep.
    "paper_q125": dict(W0=5e-3, Dn_max=-1e-4, node=5303, Q=1.25e-4),
    "paper_q135": dict(W0=5e-3, Dn_max=-1e-4, node=5303, Q=1.35e-4),
    "paper_q150": dict(W0=5e-3, Dn_max=-1e-4, node=5303, Q=1.50e-4),
    "paper_q175": dict(W0=5e-3, Dn_max=-1e-4, node=5303, Q=1.75e-4),
    # Values actually shipped in the author's first archive.
    "source5203_q100": dict(W0=1e-3, Dn_max=-1e-3, node=5203, Q=1.00e-4),
    "source5303_q100": dict(W0=1e-3, Dn_max=-1e-3, node=5303, Q=1.00e-4),
    "source5303_q125": dict(W0=1e-3, Dn_max=-1e-3, node=5303, Q=1.25e-4),
    "source5303_q150": dict(W0=1e-3, Dn_max=-1e-3, node=5303, Q=1.50e-4),
    # Intermediate aperture used to determine the hydraulic sensitivity direction.
    "mid5303_q100": dict(W0=2e-3, Dn_max=-1e-3, node=5303, Q=1.00e-4),
}


def replace_one(text: str, pattern: str, replacement: str, label: str) -> str:
    updated, count = re.subn(pattern, replacement, text, count=1, flags=re.MULTILINE)
    if count != 1:
        raise RuntimeError(f"{label}: expected one match, found {count}")
    return updated


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-dir", type=Path, required=True)
    parser.add_argument("--variant", choices=sorted(VARIANTS), required=True)
    args = parser.parse_args()
    cfg = VARIANTS[args.variant]

    params_path = args.case_dir / "defineModelParameters.m"
    params = params_path.read_text(encoding="utf-8")
    params = replace_one(
        params,
        r"params\.Dn_max\s*=\s*[-+0-9.eE]+;",
        f"params.Dn_max = {cfg['Dn_max']:.12g};",
        "Dn_max",
    )
    params = replace_one(
        params,
        r"params\.W0\s*=\s*[-+0-9.eE]+;",
        f"params.W0 = {cfg['W0']:.12g};",
        "W0",
    )
    params_path.write_text(params, encoding="utf-8")

    bc_path = args.case_dir / "defineBoundaryConditions.m"
    bc = bc_path.read_text(encoding="utf-8")
    bc = replace_one(
        bc,
        r"injectionNode\s*=\s*find\(abs\(R\(:,1\)-50\)<1e-12\s*&\s*abs\(R\(:,2\)-52\)<1e-12,\s*1\);",
        f"injectionNode = {cfg['node']};",
        "injection node",
    )
    bc = replace_one(
        bc,
        r"boundaryConditions\.concentratedSources\s*=\s*\[injectionNode,\s*[-+0-9.eE]+\];",
        f"boundaryConditions.concentratedSources = [injectionNode,{cfg['Q']:.12g}];",
        "injection rate",
    )
    bc_path.write_text(bc, encoding="utf-8")

    payload = {"variant": args.variant, **cfg}
    (args.case_dir / "VARIANT.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    (args.case_dir / "REPRODUCTION_PATCH.txt").write_text(
        f"calibration variant={args.variant}; W0={cfg['W0']:.12g}; "
        f"Dn_max={cfg['Dn_max']:.12g}; node={cfg['node']}; Q={cfg['Q']:.12g}; "
        "horizon=600h; event threshold=1e-4m/s\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, sort_keys=True))


if __name__ == "__main__":
    main()
