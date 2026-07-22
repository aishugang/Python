#!/usr/bin/env python3
"""Apply direct initial-stress states to a prepared Table 5 MATLAB case."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def replace_once(text: str, old: str, new: str, label: str) -> str:
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{label}: expected one match, found {count}")
    return text.replace(old, new)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case-dir", type=Path, required=True)
    ap.add_argument("--variant", choices=["bulk_only", "bulk_fault_projected"], required=True)
    args = ap.parse_args()
    case = args.case_dir

    params_path = case / "defineModelParameters.m"
    params = params_path.read_text(encoding="utf-8")
    params = replace_once(params, "params.Sxxr=0;", "params.Sxxr=-40e6;", "Sxxr")
    params = replace_once(params, "params.Syyr=0;", "params.Syyr=-20e6;", "Syyr")
    params = replace_once(params, "params.Sxyr=0;", "params.Sxyr=0.0;", "Sxyr")
    params_path.write_text(params, encoding="utf-8")

    main_path = case / "XFEM_MainCode.m"
    code = main_path.read_text(encoding="utf-8")
    code = replace_once(
        code,
        "params = defineModelParameters();",
        "params = defineModelParameters();\n"
        "sigma0 = [params.Sxxr, params.Sxyr; params.Sxyr, params.Syyr];\n"
        "Sn0 = n' * sigma0 * n;\n"
        "Ss0 = t' * sigma0 * n;\n"
        "fprintf('DIRECT_INITIAL_STRESS sigma_xx=%.12g sigma_yy=%.12g sigma_xy=%.12g Sn0=%.12g Ss0=%.12g\\n', params.Sxxr, params.Syyr, params.Sxyr, Sn0, Ss0);",
        "project initial stress",
    )

    if args.variant == "bulk_only":
        # Only the bulk initial-stress field is added through InsituStress.
        # Keep the author's original fault-history initialization unchanged.
        marker = "U = zeros(DOF_u+DOF_lag, 1); "
        replacement = marker + "\nfprintf('DIRECT_INITIAL_VARIANT bulk_only LM0=0 Sn_history0=%.12g Ss_history0=%.12g\\n', Sn(1,1), Ss(1,1));"
        code = replace_once(code, marker, replacement, "bulk-only marker")
    else:
        init_old = """W (:,1) = params.W0;
uf(:,1) = params.uf0;
theta(:,1) = params.theta0;
P(:,1)=params.P0;
Time(1,1)=0;
dampingFactor = 1;
Ss(:,1)=0;
Sn(:,1)=-3e+7;
Ds(:,1)=0;"""
        init_new = """uf(:,1) = params.uf0;
theta(:,1) = params.theta0;
P(:,1)=params.P0;
Time(1,1)=0;
dampingFactor = 1;
Sn(:,1)=Sn0;
Ss(:,1)=Ss0;
Dn(:,1)=0;
Ds(:,1)=0;
slip(:,1)=0;
Vs(:,1)=0;
stickflag(:,1)=0;
duf(:)=0;
W(:,1) = params.W0 + min(0,Sn(:,1))*params.Dn_max ./ (params.Kn*params.Dn_max + min(0,Sn(:,1)));"""
        code = replace_once(code, init_old, init_new, "projected fault history")

        u_old = """U = zeros(DOF_u+DOF_lag, 1); 
U_ddot = zeros(DOF_u+DOF_lag, 1); 
U_dot = zeros(DOF_u+DOF_lag, 1); """
        u_new = """U = zeros(DOF_u+DOF_lag, 1); 
U_ddot = zeros(DOF_u+DOF_lag, 1); 
U_dot = zeros(DOF_u+DOF_lag, 1); 
LM(:,1) = Sn0;
U(DOF_u+1:DOF_u+DOF_lag,1) = LM(:,1);
fprintf('DIRECT_INITIAL_VARIANT bulk_fault_projected LM0=%.12g Sn_history0=%.12g Ss_history0=%.12g W0eff=%.12g\\n', LM(1,1), Sn(1,1), Ss(1,1), W(1,1));"""
        code = replace_once(code, u_old, u_new, "initialize LM block")

    residual_marker = "%% Newton-Raphson iteration\n"
    residual_block = r'''%% Direct initial-state residual audit
initialFaultForce = FaultForce(connectivity, node_coordinates, NODES, CRACK, Pn(:,1), Ss(:,1), t, n, elementType, PHI, omega, domain_length, nex);
initialResidual = globalK * U(:,1) + [initialFaultForce; zeros(DOF_lag,1)] - globalFu;
initialResidual(dispDOFs) = 0;
initialResidualAbs = norm(initialResidual);
initialResidualScale = max(norm(globalFu_ext1) + norm(globalFu_ext3), 1);
initialResidualRatio = initialResidualAbs / initialResidualScale;
initialMaxU = max(abs(U(1:DOF_u,1)));
initialMaxLM = max(abs(LM(:,1)));
initialSummary = table(initialResidualAbs, initialResidualScale, initialResidualRatio, initialMaxU, initialMaxLM, Sn0, Ss0, ...
    'VariableNames', {'residual_abs','residual_scale','residual_ratio','max_abs_U','max_abs_LM','Sn0_Pa','Ss0_Pa'});
writetable(initialSummary, 'DirectInitialState_summary.csv');
save('DirectInitialState.mat', 'initialSummary', 'sigma0', 'Sn0', 'Ss0', 'U', 'LM', 'Sn', 'Ss', 'W', '-v7');
fprintf('DIRECT_INITIAL_RESIDUAL abs=%.12g scale=%.12g ratio=%.12g maxU=%.12g maxLM=%.12g\n', initialResidualAbs, initialResidualScale, initialResidualRatio, initialMaxU, initialMaxLM);

%% Newton-Raphson iteration
'''
    code = replace_once(code, residual_marker, residual_block, "initial residual audit")
    main_path.write_text(code, encoding="utf-8")

    payload = {
        "variant": args.variant,
        "Sxxr": -40e6,
        "Syyr": -20e6,
        "Sxyr": 0.0,
        "boundary_tractions_retained": True,
        "Q": 1e-4,
        "horizon_h": 600,
    }
    (case / "DIRECT_INITIAL_VARIANT.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    (case / "DIRECT_INITIAL_PATCH.txt").write_text(
        f"variant={args.variant}; direct bulk initial stress via InsituStress; published boundary tractions retained; Q=1e-4; 600h\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, sort_keys=True))


if __name__ == "__main__":
    main()
