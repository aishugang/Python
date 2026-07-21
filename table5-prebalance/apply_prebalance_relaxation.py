#!/usr/bin/env python3
"""Patch the prepared Table 5 case with a load-ramped drained prebalance.

The author's original coupled Newton/Newmark loop is reused. During the first
80 steps injection is disabled, pore pressure is constrained to zero, and the
published mechanical boundary load is ramped over 20 steps then held for 60
steps. At the end of the prebalance phase the equilibrated displacement/contact
state is retained, while velocity, acceleration, slip velocity, cumulative slip,
RSF state history and the injection clock are reset before Q=1e-4 is enabled.
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
    path = args.case_dir / "XFEM_MainCode.m"
    code = path.read_text(encoding="utf-8")

    code = replace_once(
        code,
        "Nt = 10000;",
        "prebalanceSteps = 80;\nprebalanceRampSteps = 20;\nNt = 10000 + prebalanceSteps;\ninjectionStartIndex = prebalanceSteps + 1;",
        "prebalance counters",
    )

    code = replace_once(
        code,
        "globalFu = [globalFu; sparse(DOF_lag, 1)];",
        "globalFu = [globalFu; sparse(DOF_lag, 1)];\nglobalFuPublished = globalFu;\nglobalFpInjection = globalFp_ext3;",
        "save published loads",
    )

    code = replace_once(
        code,
        "for nt = 1:Nt\n",
        "for nt = 1:Nt\n\n"
        "    prebalancePhase = (nt <= prebalanceSteps);\n"
        "    if prebalancePhase\n"
        "        loadFactor = min(nt / prebalanceRampSteps, 1);\n"
        "        globalFu = loadFactor * globalFuPublished;\n"
        "        globalFp_ext3 = zeros(size(globalFpInjection));\n"
        "        if nt == 1 || nt == prebalanceRampSteps || nt == prebalanceSteps\n"
        "            fprintf('PREBALANCE_PHASE step=%d loadFactor=%.12g\\n', nt, loadFactor);\n"
        "        end\n"
        "    else\n"
        "        globalFu = globalFuPublished;\n"
        "        globalFp_ext3 = globalFpInjection;\n"
        "    end\n",
        "phase switch",
    )

    adaptive_old = """    if all(Vs(:,nt+1) == 1e-9)
        dt_ev = inf;
    else
        dt_ev = min(0.1 * params.Dc ./ Vs(:,nt+1));
    end
    dt_candidate = min([max(1e-3, dt_ev), dt_c]);
    dt = min(dt_candidate, 1.2 * dt);"""
    adaptive_new = """    if prebalancePhase
        dt = dt_c;
    else
        if all(Vs(:,nt+1) == 1e-9)
            dt_ev = inf;
        else
            dt_ev = min(0.1 * params.Dc ./ Vs(:,nt+1));
        end
        dt_candidate = min([max(1e-3, dt_ev), dt_c]);
        dt = min(dt_candidate, 1.2 * dt);
    end"""
    code = replace_once(code, adaptive_old, adaptive_new, "prebalance fixed dt")

    pressure_marker = """        if ~isempty(pressureDOFs)
        Jacobian(pressureDOFs, :) = 0;
        Jacobian(:, pressureDOFs) = 0;
        diagPressIdx = sub2ind(size(Jacobian), pressureDOFs, pressureDOFs);
        Jacobian(diagPressIdx) = 1;
        Residual(pressureDOFs) = P(fixedPressureDOF(:,1), nt+1) - fixedPressures;
        end


Jt = sparse(Jacobian);"""
    pressure_repl = """        if ~isempty(pressureDOFs)
        Jacobian(pressureDOFs, :) = 0;
        Jacobian(:, pressureDOFs) = 0;
        diagPressIdx = sub2ind(size(Jacobian), pressureDOFs, pressureDOFs);
        Jacobian(diagPressIdx) = 1;
        Residual(pressureDOFs) = P(fixedPressureDOF(:,1), nt+1) - fixedPressures;
        end

        if prebalancePhase
            allPressureDOFs = (DOF_u + DOF_lag) + (1:DOF_p);
            Jacobian(allPressureDOFs, :) = 0;
            Jacobian(:, allPressureDOFs) = 0;
            diagAllP = sub2ind(size(Jacobian), allPressureDOFs, allPressureDOFs);
            Jacobian(diagAllP) = 1;
            Residual(allPressureDOFs) = P(:,nt+1);
        end

Jt = sparse(Jacobian);"""
    code = replace_once(code, pressure_marker, pressure_repl, "drained pressure constraint")

    time_old = """Time(nt+1,1) = dt + Time(nt,1);
current_time_hours = Time(end) / 3600;
if current_time_hours >= 600
        fprintf('Reached total simulation time of 600 hours. Quitting simulation.\n');
        break;
end"""
    time_new = """Time(nt+1,1) = dt + Time(nt,1);

if prebalancePhase && nt == prebalanceSteps
    prebalanceFinalNorm = Norm;
    prebalanceFinalIter = iter;
    prebalanceMaxU = max(abs(U(1:DOF_u,nt+1)));
    prebalanceMaxV = max(abs(U_dot(1:DOF_u,nt+1)));
    prebalanceFrictionRatio = max(abs(Ss(:,nt+1)) ./ max(params.uf0*abs(Sn(:,nt+1)), eps));
    preSummary = table(prebalanceFinalIter, prebalanceFinalNorm, prebalanceMaxU, prebalanceMaxV, min(Sn(:,nt+1)), max(Sn(:,nt+1)), min(Ss(:,nt+1)), max(Ss(:,nt+1)), prebalanceFrictionRatio, ...
        'VariableNames', {'iterations_last_step','last_step_norm','max_abs_U','max_abs_velocity_before_reset','min_Sn_Pa','max_Sn_Pa','min_Ss_Pa','max_Ss_Pa','max_friction_ratio'});
    writetable(preSummary, 'Prebalance_summary.csv');
    save('Prebalance_state.mat', 'preSummary', 'U', 'LM', 'Dn', 'Ds', 'Sn', 'Ss', 'W', 'P', '-v7');
    fprintf('PREBALANCE_COMPLETE norm=%.12g maxU=%.12g maxV=%.12g ratio=%.12g\\n', prebalanceFinalNorm, prebalanceMaxU, prebalanceMaxV, prebalanceFrictionRatio);

    U_dot(:,nt+1) = 0;
    U_ddot(:,nt+1) = 0;
    U_dotx(:,nt+1) = 0; U_doty(:,nt+1) = 0; U_dot_total(:,nt+1) = 0;
    U_ddotx(:,nt+1) = 0; U_ddoty(:,nt+1) = 0; U_ddot_total(:,nt+1) = 0;
    P(:,nt+1) = params.P0; Pp(:,nt+1) = params.P0;
    Time(1:nt+1,1) = 0;
    slip(:,1:nt+1) = 0;
    Vs(:,1:nt+1) = 0;
    theta(:,1:nt+1) = params.theta0;
    uf(:,1:nt+1) = params.uf0;
    stickflag(:,1:nt+1) = 0;
    duf(:) = 0;
    fprintf('PREBALANCE_DYNAMIC_STATE_RESET injection_clock_h=0 Q=%.12g\\n', globalFpInjection(find(globalFpInjection,1)));
end

current_time_hours = Time(nt+1) / 3600;
if ~prebalancePhase && current_time_hours >= 600
        fprintf('Reached total simulation time of 600 hours. Quitting simulation.\n');
        break;
end"""
    code = replace_once(code, time_old, time_new, "reset and injection clock")

    path.write_text(code, encoding="utf-8")
    (args.case_dir / "PREBALANCE_PATCH.txt").write_text(
        "80-step drained prebalance using original coupled Newton/Newmark loop; "
        "published boundary traction ramped over 20 steps; Q=0 and P=0; "
        "retain equilibrated U/LM/contact state and reset dynamics/RSF clock before injection\n",
        encoding="utf-8",
    )
    print(f"Applied load-ramped prebalance patch to {path}")


if __name__ == "__main__":
    main()
