#!/usr/bin/env python3
"""Insert a drained mechanical pre-equilibration stage before injection.

The published boundary tractions remain active. The prebalance solve keeps
P=0 and Q=0, solves the static displacement/contact equilibrium, transfers
U/LM/Sn/Ss/Dn/Ds/W into the injection initial state, and resets velocity,
acceleration, slip velocity, and simulation time before the 600 h run.
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
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-dir", type=Path, required=True)
    args = parser.parse_args()

    main_path = args.case_dir / "XFEM_MainCode.m"
    code = main_path.read_text(encoding="utf-8")

    marker = "%% Newton-Raphson iteration\n"
    block = r'''%% Drained mechanical pre-equilibration before injection
% Keep the published right/top boundary tractions active, but solve a static
% P=0 equilibrium before the injection clock starts. The converged
% displacement/contact state is retained; all dynamic state is reset.
prebalanceTolerance = 1e-9;
prebalanceMaxIter = 80;
prebalanceResidualHistory = nan(prebalanceMaxIter,1);
prebalanceU = U(:,1);
prebalanceSn = Sn(:,1);
prebalanceSs = Ss(:,1);
prebalanceDn = Dn(:,1);
prebalanceDs = Ds(:,1);
prebalanceW = W(:,1);
prebalanceLM = LM(:,1);
prebalanceStick = zeros(nInterfaceElem,1);
prebalanceDuf = zeros(nInterfaceElem,1);
prebalanceDt = dt_c; % tangent regularization only; no pseudo-time is advanced

fprintf('PREBALANCE_BEGIN mode=drained_static_reset_dynamics Q=0 P=0\n');
for preIter = 1:prebalanceMaxIter
    preKinterface = StiffnessInterface_Lagrange(params, connectivity, node_coordinates, NODES, CRACK, n, t, prebalanceSn, prebalanceDuf, prebalanceDt, prebalanceStick, elementType, PHI, omega, domain_length, nex);
    preFaultForce = FaultForce(connectivity, node_coordinates, NODES, CRACK, Pn(:,1), prebalanceSs, t, n, elementType, PHI, omega, domain_length, nex);

    preK = globalK;
    preK(1:DOF_u,1:DOF_u) = preK(1:DOF_u,1:DOF_u) + preKinterface;
    preResidual = preK * prebalanceU + [preFaultForce; zeros(DOF_lag,1)] - globalFu;

    % Essential displacement conditions are homogeneous in Example 4.
    preK(dispDOFs,:) = 0;
    preK(:,dispDOFs) = 0;
    preDiag = sub2ind(size(preK), dispDOFs, dispDOFs);
    preK(preDiag) = 1;
    preResidual(dispDOFs) = prebalanceU(dispDOFs);

    preNorm = norm(preResidual) / max(norm(globalFu), 1);
    prebalanceResidualHistory(preIter) = preNorm;
    fprintf('PREBALANCE_ITER iter=%d norm=%.12g\n', preIter, preNorm);
    if preNorm < prebalanceTolerance
        break;
    end

    preDelta = -sparse(preK) \ preResidual;
    prebalanceU = prebalanceU + preDelta;

    prebalanceLM = prebalanceU(DOF_u+1:end);
    [~, prePt, preGn, preGt, preLMelem] = Interface(params, connectivity, node_coordinates, NODES, CRACK, prebalanceU, n, t, interface_elements, prebalanceLM, elementType, omega, domain_length, nex, PHI);
    prebalanceDn = preGn;
    prebalanceDs = preGt;
    prebalanceSn = preLMelem;
    prebalanceSs = prePt;
    prebalanceW = params.W0 + min(0, prebalanceSn) .* params.Dn_max ./ (params.Kn * params.Dn_max + min(0, prebalanceSn));
end

if preNorm >= prebalanceTolerance
    error('Prebalance failed to converge: norm=%.12g after %d iterations', preNorm, preIter);
end

% Transfer the equilibrated configuration to injection t=0.
U(:,1) = prebalanceU;
LM(:,1) = prebalanceLM;
Dn(:,1) = prebalanceDn;
Ds(:,1) = prebalanceDs;
Sn(:,1) = prebalanceSn;
Ss(:,1) = prebalanceSs;
W(:,1) = prebalanceW;
Ux(:,1) = U(1:2:2*nNodes,1);
Uy(:,1) = U(2:2:2*nNodes,1);
a(:,1) = U(2*nNodes+1:2*nNodes+2*nnz(NODES(:,2)),1);

% Reset only the dynamic/injection clock state. Preserve the equilibrated
% displacement and contact state; keep initial pore pressure at zero.
U_dot(:) = 0;
U_ddot(:) = 0;
U_dotx(:) = 0; U_doty(:) = 0; U_dot_total(:) = 0;
U_ddotx(:) = 0; U_ddoty(:) = 0; U_ddot_total(:) = 0;
P(:,1) = params.P0;
Pp(:,1) = params.P0;
Time(1,1) = 0;
slip(:,1) = 0;
Vs(:,1) = 0;
theta(:,1) = params.theta0;
uf(:,1) = params.uf0;
stickflag(:,1) = 0;
duf(:) = 0;

preMu = params.uf0 + params.a*log(params.V0/params.V0) + params.b*log(params.V0*params.theta0/params.Dc);
preFrictionRatio = max(abs(Ss(:,1)) ./ max(preMu*abs(Sn(:,1)), eps));
preSummary = table(preIter, preNorm, max(abs(U(1:DOF_u,1))), min(Sn(:,1)), max(Sn(:,1)), min(Ss(:,1)), max(Ss(:,1)), preFrictionRatio, ...
    'VariableNames', {'iterations','relative_residual','max_abs_U','min_Sn_Pa','max_Sn_Pa','min_Ss_Pa','max_Ss_Pa','max_friction_ratio'});
writetable(preSummary, 'Prebalance_summary.csv');
save('Prebalance_state.mat', 'preSummary', 'prebalanceResidualHistory', 'U', 'LM', 'Dn', 'Ds', 'Sn', 'Ss', 'W', 'P', '-v7');
fprintf('PREBALANCE_COMPLETE iter=%d norm=%.12g maxU=%.12g Sn=[%.12g,%.12g] Ss=[%.12g,%.12g] ratio=%.12g\n', ...
    preIter, preNorm, max(abs(U(1:DOF_u,1))), min(Sn(:,1)), max(Sn(:,1)), min(Ss(:,1)), max(Ss(:,1)), preFrictionRatio);

%% Newton-Raphson iteration
'''
    code = replace_once(code, marker, block, "prebalance insertion")
    main_path.write_text(code, encoding="utf-8")

    (args.case_dir / "PREBALANCE_PATCH.txt").write_text(
        "drained static mechanical pre-equilibration with published boundary tractions; "
        "P=0 and injection disabled during equilibrium; retain U/LM/contact state; "
        "reset velocities, accelerations, slip velocity and Time before Q=1e-4 injection\n",
        encoding="utf-8",
    )
    print(f"Applied prebalance-reset patch to {main_path}")


if __name__ == "__main__":
    main()
