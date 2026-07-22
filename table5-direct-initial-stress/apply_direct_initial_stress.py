#!/usr/bin/env python3
"""Apply a direct initial in-situ stress/contact state to the prepared Table 5 case.

No prebalance stage is used. The initial continuum stress is assembled through
InsituStress, consistent right/top tractions are retained, and the 45-degree
fault is initialized with the projected normal/tangential traction. Physical
U, velocity and acceleration remain zero; the LM block carries the initial
normal contact pressure.
"""
from __future__ import annotations

import argparse
from pathlib import Path


def replace_once(text: str, old: str, new: str, label: str) -> str:
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{label}: expected one match, found {count}")
    return text.replace(old, new)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case-dir", type=Path, required=True)
    args = ap.parse_args()

    params_path = args.case_dir / "defineModelParameters.m"
    params = params_path.read_text(encoding="utf-8")
    params = replace_once(params, "params.Sxxr=0;", "params.Sxxr=-40e6;", "Sxxr")
    params = replace_once(params, "params.Syyr=0;", "params.Syyr=-20e6;", "Syyr")
    params = replace_once(params, "params.Sxyr=0;", "params.Sxyr=0.0;", "Sxyr")
    params_path.write_text(params, encoding="utf-8")

    main_path = args.case_dir / "XFEM_MainCode.m"
    code = main_path.read_text(encoding="utf-8")

    old_initial = """W (:,1) = params.W0;
uf(:,1) = params.uf0;
theta(:,1) = params.theta0;
P(:,1)=params.P0;
Time(1,1)=0;
dampingFactor = 1;
Ss(:,1)=0;
Sn(:,1)=-3e+7;
Ds(:,1)=0;"""
    new_initial = """% Direct initial stress projected onto the fault.
sigma0 = [params.Sxxr, params.Sxyr; params.Sxyr, params.Syyr];
Sn0 = n' * sigma0 * n;
Ss0 = t' * sigma0 * n;

uf(:,1) = params.uf0;
theta(:,1) = params.theta0;
P(:,1) = params.P0;
Pp(:,1) = params.P0;
Time(1,1) = 0;
dampingFactor = 1;
Sn(:,1) = Sn0;
Ss(:,1) = Ss0;
Dn(:,1) = 0;
Ds(:,1) = 0;
slip(:,1) = 0;
Vs(:,1) = params.V0;
stickflag(:,1) = 0;
duf(:) = 0;
W(:,1) = params.W0 + min(0,Sn(:,1))*params.Dn_max ./ ...
    (params.Kn*params.Dn_max + min(0,Sn(:,1)));
fprintf('DIRECT_INITIAL_FAULT_TRACTION Sn0=%.12g Ss0=%.12g\\n', Sn0, Ss0);"""
    code = replace_once(code, old_initial, new_initial, "initial state")

    old_u = """U = zeros(DOF_u+DOF_lag, 1); 
U_ddot = zeros(DOF_u+DOF_lag, 1); 
U_dot = zeros(DOF_u+DOF_lag, 1); 
ResidualHistory = cell(Nt+1,1);"""
    new_u = """U = zeros(DOF_u+DOF_lag, 1);
U_ddot = zeros(DOF_u+DOF_lag, 1);
U_dot = zeros(DOF_u+DOF_lag, 1);

% LM is an independent unknown and must carry the initial normal traction.
LM(:,1) = Sn0;
U(1:DOF_u,1) = 0;
U(DOF_u+1:DOF_u+DOF_lag,1) = LM(:,1);
U_dot(:,1) = 0;
U_ddot(:,1) = 0;
ResidualHistory = cell(Nt+1,1);"""
    code = replace_once(code, old_u, new_u, "LM initial state")

    marker = "%% Newton-Raphson iteration\n"
    diagnostic = r'''%% Direct initial-state equilibrium diagnostic
initialKinterface = StiffnessInterface_Lagrange(params, connectivity, node_coordinates, NODES, CRACK, n, t, Sn(:,1), duf, dt, stickflag(:,1), elementType, PHI, omega, domain_length, nex);
initialFaultForce = FaultForce(connectivity, node_coordinates, NODES, CRACK, Pn(:,1), Ss(:,1), t, n, elementType, PHI, omega, domain_length, nex);
initialK = globalK;
initialK(1:DOF_u,1:DOF_u) = initialK(1:DOF_u,1:DOF_u) + initialKinterface;
initialResidual = initialK*U(:,1) + [initialFaultForce;zeros(DOF_lag,1)] - globalFu;
initialResidual(dispDOFs) = 0;
initialResidualRatio = norm(initialResidual)/max(norm(globalFu),1);
initialFrictionRatio = max(abs(Ss(:,1))./max(params.uf0*abs(Sn(:,1)),eps));
initialSummary = table(Sn0,Ss0,initialResidualRatio,initialFrictionRatio,max(abs(U(1:DOF_u,1))), ...
    'VariableNames',{'Sn0_Pa','Ss0_Pa','relative_mechanical_residual','max_friction_ratio','max_abs_physical_U'});
writetable(initialSummary,'Direct_initial_stress_summary.csv');
save('Direct_initial_stress_state.mat','initialSummary','U','LM','Sn','Ss','Dn','Ds','W','P','-v7');
fprintf('DIRECT_INITIAL_STATE residual=%.12g friction_ratio=%.12g maxPhysicalU=%.12g\n', initialResidualRatio, initialFrictionRatio, max(abs(U(1:DOF_u,1))));

%% Newton-Raphson iteration
'''
    code = replace_once(code, marker, diagnostic, "initial residual diagnostic")
    main_path.write_text(code, encoding="utf-8")

    (args.case_dir / "DIRECT_INITIAL_STRESS_PATCH.txt").write_text(
        "Direct t=0 initial state: sigma_xx=-40 MPa, sigma_yy=-20 MPa, sigma_xy=0; "
        "retain matching right/top tractions; project sigma0 to fault Sn0/Ss0; "
        "initialize LM with Sn0; U/Udot/Uddot/P/slip remain zero; no prebalance.\n",
        encoding="utf-8",
    )
    print(f"Applied direct initial stress state to {main_path}")


if __name__ == "__main__":
    main()
