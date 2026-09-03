%VALIDATE_CELL_RANGE Rigorous, fail-closed Omega_mid coverage certificate.
%
% This compatibility entry point remains a script, so existing invocations
% such as run('tools/validate_cell_range.m') continue to work.  The proof
% implementation imports CSV endpoints as text and converts their decimal
% strings directly to INTLAB intervals; it uses no geometric tolerance.

coverage_project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(coverage_project_root);
addpath(fullfile(coverage_project_root,'scripts_run'));
addpath(fullfile(coverage_project_root,'tools'));

coverage_input_file = fullfile( ...
    coverage_project_root,'inputs','cell_def.csv');
coverage_manifest_file = fullfile( ...
    coverage_project_root,'results', ...
    'omega_mid_geometric_coverage_manifest.json');

global INTERVAL_MODE
if isempty(INTERVAL_MODE) || ~logical(INTERVAL_MODE)
    omega_up_all_prepare_worker(coverage_project_root,'interval');
end

coverage_certificate = validate_omega_mid_coverage( ...
    coverage_input_file,'manifest_file',coverage_manifest_file);
if ~coverage_certificate.complete
    if isempty(coverage_certificate.witness)
        coverage_failure_kind = 'unspecified';
    else
        coverage_failure_kind = coverage_certificate.witness.kind;
    end
    error('validate_cell_range:CoverageFailed', ...
        'Omega_mid coverage failed (%s).',coverage_failure_kind);
end

fprintf('OK: cell_def.csv rigorously covers Omega_mid.\n');
fprintf('Coverage proof payload SHA-256: %s\n', ...
    coverage_certificate.hashes.proof_payload);
fprintf('Coverage manifest: %s\n',coverage_manifest_file);
