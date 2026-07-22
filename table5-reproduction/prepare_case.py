#!/usr/bin/env python3
"""Prepare a headless MATLAB R2022b reproduction of Sabah et al. Table 5.

Starting point: the author's first public MATLAB archive (commit 8f7b137).
This baseline follows the paper literally: 600 h, W0=5 mm,
Dn_max=-0.1 mm, injection at the mesh node (50, 52), and seismic-event
threshold 1e-4 m/s.  Post-processing uses event-increment slip.
"""
from __future__ import annotations

import argparse
import shutil
import zipfile
from pathlib import Path


def replace_exact(text: str, old: str, new: str, label: str) -> str:
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{label}: expected one match, found {count}")
    return text.replace(old, new)


def write_headless_mesh(dst: Path) -> None:
    (dst / "MeshGeneration_2D.m").write_text(r'''function [connectivity_global, node_coordinates_global, nodes, elements] = MeshGeneration_2D(domain_length, domain_width, nex, ney, crackStart, crackEnd)
% Exact structured-mesh numbering used by the author, without prompts/plots.
x = linspace(0, domain_length, nex + 1);
y = linspace(0, domain_width, ney + 1);
[X, Y] = meshgrid(x, y);
X = X'; Y = Y';
node_coordinates_global = [X(:), Y(:)];
nodes = size(node_coordinates_global, 1);
elements = nex * ney;
connectivity_global = zeros(elements, 4);
el_id = 1;
for j = 1:ney
    for i = 1:nex
        n1 = (j - 1) * (nex + 1) + i;
        n2 = n1 + 1;
        n3 = n2 + nex + 1;
        n4 = n1 + nex + 1;
        connectivity_global(el_id, :) = [n1, n2, n3, n4];
        el_id = el_id + 1;
    end
end
if isempty(crackStart) || isempty(crackEnd)
    error('Crack endpoints must be supplied.');
end
end
''', encoding="utf-8")


def write_metrics(dst: Path) -> None:
    (dst / "Table5Metrics.m").write_text(r'''function metrics = Table5Metrics(interface_elements, interface_nodes, Vs, slip, Ss, Sn, params, Time)
threshold = 1e-4;                  % paper Section 4.4 [m/s]
faultWidth = 1.0;                  % 2-D unit thickness [m]
numElems = size(interface_elements, 1);
elementLength = zeros(numElems, 1);
for e = 1:numElems
    p1 = interface_nodes(interface_elements(e,1), :);
    p2 = interface_nodes(interface_elements(e,2), :);
    elementLength(e) = norm(p2 - p1);
end
faultLength = sum(elementLength);
active = any(abs(Vs) > threshold, 1);
transitions = diff([false, active, false]);
starts = find(transitions == 1);
ends = find(transitions == -1) - 1;
numEvents = min(numel(starts), numel(ends));
metrics = zeros(numEvents, 7);
fid = fopen('Table5_reproduced.csv', 'w');
assert(fid >= 0, 'Cannot open Table5_reproduced.csv');
fprintf(fid, 'event,onset_h,average_slip_m,seismic_moment_Nm,moment_magnitude,static_stress_drop_Pa,max_slip_velocity_mps,rupture_duration_s\n');
for k = 1:numEvents
    i0 = starts(k); i1 = ends(k); iprev = max(1, i0 - 1);
    eventSlip = abs(slip(:,i1) - slip(:,iprev));
    averageSlip = sum(eventSlip .* elementLength) / faultLength;
    seismicMoment = params.G * faultWidth * sum(eventSlip .* elementLength);
    momentMagnitude = (2/3) * (log10(max(seismicMoment, eps)) - 9.05);
    tauBefore = abs(Ss(:,iprev)); tauAfter = abs(Ss(:,i1));
    stressDrop = sum((tauBefore - tauAfter) .* elementLength) / faultLength;
    onsetHours = Time(i0) / 3600;
    durationSeconds = max(0, Time(i1) - Time(i0));
    maxVelocity = max(max(abs(Vs(:,i0:i1))));
    metrics(k,:) = [onsetHours, averageSlip, seismicMoment, momentMagnitude, stressDrop, maxVelocity, durationSeconds];
    fprintf(fid, '%d,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n', k, metrics(k,:));
    fprintf('TABLE5_EVENT %d onset_h=%.9g avg_slip=%.9g M0=%.9g Mw=%.9g stress_drop=%.9g maxV=%.9g duration_s=%.9g\n', k, metrics(k,:));
end
fclose(fid);
save('Table5_results.mat', 'metrics', 'starts', 'ends', 'threshold', 'faultLength', 'Time', 'Vs', 'slip', 'Ss', 'Sn', '-v7');
end
''', encoding="utf-8")


def write_runner(dst: Path) -> None:
    (dst / "run_table5_ci.m").write_text(r'''function run_table5_ci
% Preserve a complete diagnostic log even when the author code fails.
diary('run.log');
fprintf('MATLAB_VERSION %s\n', version);
fprintf('CASE_DIRECTORY %s\n', pwd);
try
    XFEM_MainCode;
    fprintf('MATLAB_CASE_SUCCESS\n');
    diary off;
catch ME
    fprintf(2, 'MATLAB_CASE_FAILURE %s\n', ME.message);
    fid = fopen('run_error.txt', 'w');
    if fid >= 0
        fprintf(fid, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
        fclose(fid);
    end
    diary off;
    rethrow(ME);
end
end
''', encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--author-dir', type=Path, default=Path('author'))
    ap.add_argument('--output-dir', type=Path, default=Path('work/Table5_case'))
    args = ap.parse_args()

    archives = list(args.author_dir.rglob('XFEM_Induced Seismicity.zip'))
    if len(archives) != 1:
        raise RuntimeError(f'Expected one author ZIP, found {archives}')
    extract_root = args.output_dir.parent / '_extract'
    shutil.rmtree(extract_root, ignore_errors=True)
    shutil.rmtree(args.output_dir, ignore_errors=True)
    extract_root.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archives[0]) as zf:
        zf.extractall(extract_root)
    roots = [p.parent for p in extract_root.rglob('XFEM_MainCode.m')]
    if len(roots) != 1:
        raise RuntimeError(f'Expected one MATLAB source root, found {roots}')
    shutil.copytree(roots[0], args.output_dir)

    params_path = args.output_dir / 'defineModelParameters.m'
    params = params_path.read_text(encoding='utf-8-sig')
    params = replace_exact(params, 'params.Dn_max = -1e-3;', 'params.Dn_max = -1e-4;', 'paper Dn_max')
    params = replace_exact(params, 'params.W0 = 1e-3;', 'params.W0 = 5e-3;', 'paper W0')
    params_path.write_text(params, encoding='utf-8')

    bc_path = args.output_dir / 'defineBoundaryConditions.m'
    bc = bc_path.read_text(encoding='utf-8-sig')
    bc = replace_exact(
        bc,
        'boundaryConditions.concentratedSources = [5203,1e-4];',
        "injectionNode = find(abs(R(:,1)-50)<1e-12 & abs(R(:,2)-52)<1e-12, 1);\n"
        "    assert(~isempty(injectionNode), 'Paper injection node (50,52) not found');\n"
        "    boundaryConditions.concentratedSources = [injectionNode,1e-4];\n"
        "    fprintf('INJECTION_NODE id=%d x=%.12g y=%.12g\\n', injectionNode, R(injectionNode,1), R(injectionNode,2));",
        'paper injection point')
    bc_path.write_text(bc, encoding='utf-8')

    main_path = args.output_dir / 'XFEM_MainCode.m'
    code = main_path.read_text(encoding='utf-8-sig')
    code = code.replace('close all;', "% close all; % headless")
    code = replace_exact(code, 'current_time_hours >= 200', 'current_time_hours >= 600', '600 h horizon')
    code = replace_exact(code, 'Reached total simulation time of 200 hours.', 'Reached total simulation time of 600 hours.', '600 h message')
    progress_old = 'Time(nt+1,1) = dt + Time(nt,1);'
    progress_new = r'''Time(nt+1,1) = dt + Time(nt,1);
if mod(nt, 25) == 0
    fprintf('TABLE5_PROGRESS step=%d time_h=%.12g dt_s=%.12g maxV=%.12g\n', nt, Time(nt+1)/3600, dt, max(abs(Vs(:,nt+1))));
    checkpointTime = Time(1:nt+1); checkpointVs = Vs(:,1:nt+1); checkpointSlip = slip(:,1:nt+1);
    save('Table5_checkpoint.mat', 'checkpointTime', 'checkpointVs', 'checkpointSlip', '-v7');
end'''
    code = replace_exact(code, progress_old, progress_new, 'progress checkpoint')

    marker = 'Qinj = boundaryConditions.concentratedSources(1,2);'
    trim = r'''usedSteps = nt + 1;
Time = Time(1:usedSteps,:); Vs = Vs(:,1:usedSteps); slip = slip(:,1:usedSteps);
Ss = Ss(:,1:usedSteps); Sn = Sn(:,1:usedSteps); W = W(:,1:usedSteps);
Qinj = boundaryConditions.concentratedSources(1,2);'''
    code = replace_exact(code, marker, trim, 'trim histories')
    old_tail = r'''%% Calculate seismicity parameters

[AS, M0, Mw, CFS, SSD, ST] = SeismicityParameters(interface_elements, interface_nodes, Vs, slip, Ss, Sn, uf, params);

%% Calculate stresses 

[Sxx,Sxy,Syy,Svm,S1,S2,theta_p] = elemStress2(params, U, connectivity, node_coordinates, PSI, NODES, CRACK, omega, domain_length, nex, elementType, Nt);

%% Post-Processing and Visualization
postProcess(node_coordinates, U(1:2*nNodes,end), Pp(:,end), connectivity, Sxx(:,end), Syy(:,end), Sxy(:,end), elementType);
'''
    new_tail = r'''%% Table 5 event metrics
metrics = Table5Metrics(interface_elements, interface_nodes, Vs, slip, Ss, Sn, params, Time);
save('Table5_state_compact.mat', 'metrics', 'Time', 'Vs', 'slip', 'Ss', 'Sn', 'W', 'params', 'interface_elements', 'interface_nodes', 'node_coordinates', 'boundaryConditions', '-v7');
fprintf('TABLE5_RUN_COMPLETE steps=%d final_time_h=%.12g events=%d\n', usedSteps, Time(end)/3600, size(metrics,1));
'''
    code = replace_exact(code, old_tail, new_tail, 'replace graphical post-processing')
    main_path.write_text(code, encoding='utf-8')

    write_headless_mesh(args.output_dir)
    write_metrics(args.output_dir)
    write_runner(args.output_dir)
    (args.output_dir / 'REPRODUCTION_PATCH.txt').write_text(
        'paper-literal: W0=5e-3; Dn_max=-1e-4; injection=(50,52); horizon=600h; threshold=1e-4m/s\n',
        encoding='utf-8')
    print(f'Prepared {args.output_dir}')


if __name__ == '__main__':
    main()
