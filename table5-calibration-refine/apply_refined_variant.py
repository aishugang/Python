#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

VARIANTS = {
    "paper_q1100": dict(W0=5e-3, Dn_max=-1e-4, node=5303, Q=1.100e-4),
    "paper_q1125": dict(W0=5e-3, Dn_max=-1e-4, node=5303, Q=1.125e-4),
    "paper_q1150": dict(W0=5e-3, Dn_max=-1e-4, node=5303, Q=1.150e-4),
    "paper_q1175": dict(W0=5e-3, Dn_max=-1e-4, node=5303, Q=1.175e-4),
    "paper_q1200": dict(W0=5e-3, Dn_max=-1e-4, node=5303, Q=1.200e-4),
    "source_q1150": dict(W0=1e-3, Dn_max=-1e-3, node=5303, Q=1.150e-4),
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
    params = replace_one(params, r"params\.Dn_max\s*=\s*[-+0-9.eE]+;", f"params.Dn_max = {cfg['Dn_max']:.12g};", "Dn_max")
    params = replace_one(params, r"params\.W0\s*=\s*[-+0-9.eE]+;", f"params.W0 = {cfg['W0']:.12g};", "W0")
    params_path.write_text(params, encoding="utf-8")

    bc_path = args.case_dir / "defineBoundaryConditions.m"
    bc = bc_path.read_text(encoding="utf-8")
    bc = replace_one(bc, r"injectionNode\s*=\s*find\(abs\(R\(:,1\)-50\)<1e-12\s*&\s*abs\(R\(:,2\)-52\)<1e-12,\s*1\);", f"injectionNode = {cfg['node']};", "injection node")
    bc = replace_one(bc, r"boundaryConditions\.concentratedSources\s*=\s*\[injectionNode,\s*[-+0-9.eE]+\];", f"boundaryConditions.concentratedSources = [injectionNode,{cfg['Q']:.12g}];", "injection rate")
    bc_path.write_text(bc, encoding="utf-8")

    payload = {"variant": args.variant, **cfg}
    (args.case_dir / "VARIANT.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, sort_keys=True))


if __name__ == "__main__":
    main()
