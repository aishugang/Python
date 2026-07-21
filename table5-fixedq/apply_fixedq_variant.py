#!/usr/bin/env python3
"""Apply evidence-based fixed-injection-rate variants for Sabah Table 5.

All variants keep the published injection rate Q=1e-4 m^3/(s.m) and the
published well node (50,52).  The sweep isolates the inconsistent initial
stress initialization found between the first public archive and later
LM/v2 archives.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

# For a 45-degree fault under sigma_x=-40 MPa and sigma_y=-20 MPa,
# the resolved tractions are sigma_n=-30 MPa and |tau|=10 MPa.
# The sign follows the later author archive, where (-60,-20) MPa is seeded
# as Sn=-40 MPa and Ss=-20 MPa.
VARIANTS = {
    "paper_boundary_baseline": dict(
        mode="boundary", W0=5e-3, Dn=-1e-4, Sn=-30e6, Ss=0.0, Ds=0.0, theta_factor=1.0),
    "paper_boundary_ss10": dict(
        mode="boundary", W0=5e-3, Dn=-1e-4, Sn=-30e6, Ss=-10e6, Ds=0.0, theta_factor=1.0),
    "paper_boundary_ss10_ds05": dict(
        mode="boundary", W0=5e-3, Dn=-1e-4, Sn=-30e6, Ss=-10e6, Ds=-5.02413e-4, theta_factor=1.0),
    "paper_insitu_ss10": dict(
        mode="insitu", W0=5e-3, Dn=-1e-4, Sn=-30e6, Ss=-10e6, Ds=0.0, theta_factor=1.0),
    "paper_insitu_ss10_ds05": dict(
        mode="insitu", W0=5e-3, Dn=-1e-4, Sn=-30e6, Ss=-10e6, Ds=-5.02413e-4, theta_factor=1.0),
    "source_aperture_insitu_ss10_ds05": dict(
        mode="insitu", W0=1e-3, Dn=-1e-3, Sn=-30e6, Ss=-10e6, Ds=-5.02413e-4, theta_factor=1.0),
    "paper_insitu_ss8_ds04": dict(
        mode="insitu", W0=5e-3, Dn=-1e-4, Sn=-30e6, Ss=-8e6, Ds=-4.01930e-4, theta_factor=1.0),
    "paper_insitu_ss12_ds06": dict(
        mode="insitu", W0=5e-3, Dn=-1e-4, Sn=-30e6, Ss=-12e6, Ds=-6.02895e-4, theta_factor=1.0),
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
    for key, value in (("W0", cfg["W0"]), ("Dn_max", cfg["Dn"])):
        params = replace_one(
            params,
            rf"params\.{key}\s*=\s*[-+0-9.eE]+;",
            f"params.{key} = {value:.12g};",
            key,
        )

    if cfg["mode"] == "insitu":
        sxx, syy = -40e6, -20e6
    else:
        sxx, syy = 0.0, 0.0
    params = replace_one(params, r"params\.Sxxr\s*=\s*[-+0-9.eE]+;", f"params.Sxxr = {sxx:.12g};", "Sxxr")
    params = replace_one(params, r"params\.Syyr\s*=\s*[-+0-9.eE]+;", f"params.Syyr = {syy:.12g};", "Syyr")
    params = replace_one(params, r"params\.Sxyr\s*=\s*[-+0-9.eE]+;", "params.Sxyr = 0;", "Sxyr")
    params = replace_one(
        params,
        r"params\.theta0\s*=\s*params\.Dc\s*/\s*params\.V0;",
        f"params.theta0 = {cfg['theta_factor']:.12g} * params.Dc / params.V0;",
        "theta0",
    )
    params_path.write_text(params, encoding="utf-8")

    bc_path = args.case_dir / "defineBoundaryConditions.m"
    bc = bc_path.read_text(encoding="utf-8")
    if cfg["mode"] == "insitu":
        right, top = "[]", "[]"
    else:
        right, top = "[-40e+6,0]", "[0,-20e+6]"
    bc = replace_one(bc, r"tractionConfig\.rightEdge\s*=\s*[^;]+;", f"tractionConfig.rightEdge = {right};", "right traction")
    bc = replace_one(bc, r"tractionConfig\.topEdge\s*=\s*[^;]+;", f"tractionConfig.topEdge = {top};", "top traction")
    bc = replace_one(
        bc,
        r"boundaryConditions\.concentratedSources\s*=\s*\[injectionNode,\s*[-+0-9.eE]+\];",
        "boundaryConditions.concentratedSources = [injectionNode,1e-4];",
        "fixed published Q",
    )
    bc_path.write_text(bc, encoding="utf-8")

    main_path = args.case_dir / "XFEM_MainCode.m"
    code = main_path.read_text(encoding="utf-8")
    code = replace_one(code, r"Ss\(:,1\)\s*=\s*[-+0-9.eE]+;", f"Ss(:,1) = {cfg['Ss']:.12g};", "initial Ss")
    code = replace_one(code, r"Sn\(:,1\)\s*=\s*[-+0-9.eE]+;", f"Sn(:,1) = {cfg['Sn']:.12g};", "initial Sn")
    code = replace_one(code, r"Ds\(:,1\)\s*=\s*[-+0-9.eE]+;", f"Ds(:,1) = {cfg['Ds']:.12g};", "initial Ds")
    main_path.write_text(code, encoding="utf-8")

    payload = {"variant": args.variant, "Q": 1e-4, "node": 5303, **cfg}
    (args.case_dir / "FIXEDQ_VARIANT.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    (args.case_dir / "REPRODUCTION_PATCH.txt").write_text(
        "fixed published Q=1e-4 and node=(50,52); " + json.dumps(payload, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, sort_keys=True))


if __name__ == "__main__":
    main()
