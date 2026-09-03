function test_omega_mid_unified_driver()
%TEST_OMEGA_MID_UNIFIED_DRIVER  Strict LG/CR/subdivision regression.
%
% Cells 138748 and 152387 are the two rows whose historical committed CSV
% values were +Inf.  This test requires finite strict LG certificates for
% both, a direct CR certificate for cell 11731, and an adaptive subdivision
% certificate for the formerly non-closing coarse cell 1.  LG escalation
% for a CR-input row is permitted only at terminal geometric leaves.  Every
% evaluated node reuses one spectral computation for J1/J2, and a separate
% deliberately under-resolved run checks that an uncertified terminal leaf
% fails closed.

if ~strcmp(getenv('RUN_INTLAB_SMOKE'),'1')
    fprintf(['Omega_mid unified interval test skipped; set ', ...
        'RUN_INTLAB_SMOKE=1 and INTLAB_ROOT to enable it.\n']);
    return;
end

project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(project_root,'scripts_run'));
output_dir = tempname;
mkdir(output_dir);
cleanup = onCleanup(@() local_remove_tree(output_dir));
target_ids = [1,11731,138748,152387];
max_subdivision_depth = 3;
max_fem_refinement_level = 2;
minimum_cr_mesh_size = 0.005;
minimum_lg_mesh_size = 0.025;

manifest = run_omega_mid_unified_parallel( ...
    'mode','interval', ...
    'workers',0, ...
    'input_file',fullfile(project_root,'inputs','cell_def.csv'), ...
    'output_dir',output_dir, ...
    'run_name','targeted_regression', ...
    'chunk_size',2, ...
    'checkpoint_every',1, ...
    'max_subdivision_depth',max_subdivision_depth, ...
    'max_fem_refinement_level',max_fem_refinement_level, ...
    'minimum_cr_mesh_size',minimum_cr_mesh_size, ...
    'minimum_lg_mesh_size',minimum_lg_mesh_size, ...
    'cell_ids',target_ids, ...
    'resume',true, ...
    'retry_failed',true);

assert(manifest.summary.exact_selected_id_coverage);
assert(~manifest.summary.exact_full_input_id_coverage);
assert(~manifest.complete_certificate);
assert(manifest.summary.geometric_coverage_certified);
assert(manifest.geometric_coverage.complete);
assert(manifest.geometric_coverage.rigorous);
assert(strcmp(manifest.geometric_coverage.status,'certified'));
assert(manifest.geometric_coverage.proof.decimal_topology_exact);
assert(~manifest.geometric_coverage.proof.arbitrary_tolerance_used);
assert(manifest.summary.num_selected_cells == 4);
assert(manifest.summary.num_strict_records == 4);
assert(manifest.summary.num_J1_certified == 4);
assert(manifest.summary.num_J2_certified == 4);
assert(manifest.summary.num_errors == 0);
assert(manifest.summary.num_nonfinite_J1 == 0);
assert(manifest.summary.num_nonfinite_J2 == 0);
assert(manifest.summary.method_counts.CR == 1);
assert(manifest.summary.method_counts.LG == 2);
assert(manifest.summary.method_counts.SUBDIVISION == 1);
assert(manifest.summary.num_subdivided_parent_cells == 1);
assert(manifest.summary.num_fem_refined_parent_cells == 1);
assert(manifest.input.cell_count == 188623);
assert(strcmp(manifest.input.sha256,ver10_file_sha256( ...
    fullfile(project_root,'inputs','cell_def.csv'))));
assert(manifest.input.validation.ids_equal_row_numbers);
assert(numel(manifest.fingerprint) == 64);
assert(strcmp(manifest.config.schema, ...
    'lowerboundsineq.omega_mid_unified.config.v3'));
assert(isfield(manifest.config,'runtime'));
assert(isfield(manifest.config,'runtime_sha256'));
assert(numel(manifest.config.runtime_sha256) == 64);
assert(isequaln(manifest.config.runtime,manifest.runtime));
assert(numel(manifest.config.source_normalized_porcelain_sha256) == 64);
assert(numel(manifest.config.source_dirty_content_sha256) == 64);
assert(numel(manifest.config.source_state_sha256) == 64);
assert(strcmp(manifest.config.source_state_sha256, ...
    manifest.source.source_state_sha256));
assert(manifest.source.source_state_digest_captured_before_run);
assert(manifest.source.source_state_digest_rechecked_after_run);
assert(strcmp(manifest.source.source_state_sha256, ...
    manifest.source.source_state_sha256_after_run));
assert(strcmp(manifest.execution.fem_refinement_scheme, ...
    'dyadic-CR-then-terminal-LG-v1'));
assert(strcmp(manifest.execution.LG_escalation_policy, ...
    'LG-input-or-terminal-depth-after-CR-refinement-v1'));
assert(strcmp(manifest.execution.leaf_acceptance_policy, ...
    'strict-CR-or-LG-J1-and-J2-v1'));
assert(strcmp(manifest.execution.unverified_leaf_policy, ...
    'fail-closed-parent-invalid-v1'));

run_dir = fullfile(output_dir,'targeted_regression');
j1_path = fullfile(run_dir,'J1_OmegaMid.csv');
j2_path = fullfile(run_dir,'J2_OmegaMid.csv');
json_path = fullfile(run_dir,'omega_mid_unified_manifest.json');
coverage_json_path = fullfile( ...
    run_dir,'omega_mid_geometric_coverage_manifest.json');
mat_path = fullfile(run_dir,'omega_mid_unified_results.mat');
config_path = fullfile(run_dir,'config.mat');
J1 = readtable(j1_path,'TextType','string');
J2 = readtable(j2_path,'TextType','string');
saved_results = load(mat_path,'records','config');
saved_records = saved_results.records;
assert(strcmp(saved_results.config.LG_escalation_policy, ...
    'LG-input-or-terminal-depth-after-CR-refinement-v1'));
saved_config = load(config_path, ...
    'config','fingerprint','fingerprint_payload');
assert(isequaln(saved_config.config,manifest.config));
assert(strcmp(saved_config.fingerprint,manifest.fingerprint));
assert(contains(saved_config.fingerprint_payload, ...
    ['runtime=',manifest.config.runtime_sha256]));
for q = 1:numel(saved_records)
    accepted_leaves = saved_records{q}.leaves;
    assert(~isempty(accepted_leaves));
    assert(all([accepted_leaves.strict_checks_passed]));
    assert(all(ismember({accepted_leaves.method},{'CR','LG'})));
end
assert(height(J1) == 4 && height(J2) == 4);
assert(isequal(sort(double(J1.cell_id)),sort(target_ids(:))));
assert(isequal(sort(double(J2.cell_id)),sort(target_ids(:))));
assert(all(double(J1.verified) == 1));
assert(all(double(J2.verified) == 1));
assert(all(strcmpi(strtrim(string(J1.status)),'ok')));
assert(all(strcmpi(strtrim(string(J2.status)),'ok')));
assert(ismember('note',J1.Properties.VariableNames));
assert(ismember('run_timestamp',J1.Properties.VariableNames));
assert(all(isfinite(double(J1.J_lower)) & double(J1.J_lower) > 0));
assert(all(isfinite(double(J2.J_lower)) & double(J2.J_lower) > 0));
assert(all(isfinite(double(J1.lambda1_lower)) ...
    & double(J1.lambda1_lower) > 0));
assert(all(double(J1.shared_eigenvalue_computation) == 1));
assert(all(double(J1.eigenvalue_attempts) >= 1));
assert(all(strcmp(string(J1.LG_escalation_policy), ...
    "LG-input-or-terminal-depth-after-CR-refinement-v1")));
assert(all(strcmp(string(J1.fem_refinement_scheme), ...
    "dyadic-CR-then-terminal-LG-v1")));
assert(isequaln(J1.method,J2.method));
assert(isequaln(J1.lambda1_lower,J2.lambda1_lower));
assert(isequaln(J1.leaf_count,J2.leaf_count));
assert(isequaln(J1.spectral_solve_count,J2.spectral_solve_count));
assert(isequaln(J1.LG_shift_method,J2.LG_shift_method));
assert(isequaln(J1.LG_fallback_used,J2.LG_fallback_used));
assert(isequaln(J1.LG_fallback_reason,J2.LG_fallback_reason));
assert(all(isfinite(double(J1.eigenvalue_total_seconds)) ...
    & double(J1.eigenvalue_total_seconds) >= 0));
assert(all(isfinite(double(J1.LG_CR_shift_seconds)) ...
    & double(J1.LG_CR_shift_seconds) >= 0));
evaluated_nodes = double(J1.evaluated_node_count);
spectral_solves = double(J1.spectral_solve_count);
assert(all(isfinite(evaluated_nodes) ...
    & evaluated_nodes >= 1 ...
    & evaluated_nodes == floor(evaluated_nodes)));
assert(all(isfinite(spectral_solves) ...
    & spectral_solves >= evaluated_nodes ...
    & spectral_solves == floor(spectral_solves)));

k_subdivision = find(double(J1.cell_id) == 1,1);
assert(strcmp(char(J1.method(k_subdivision)),'SUBDIVISION'));
assert(double(J1.subdivision_used(k_subdivision)) == 1);
assert(double(J1.leaf_count(k_subdivision)) > 1);
assert(double(J1.max_leaf_depth(k_subdivision)) >= 1);
assert(double(J1.max_leaf_depth(k_subdivision)) ...
    == max_subdivision_depth);
assert(double(J1.leaf_count(k_subdivision)) == 64);
assert(double(J1.evaluated_node_count(k_subdivision)) ...
    == (4*double(J1.leaf_count(k_subdivision))-1)/3);
assert(double(J1.evaluated_node_count(k_subdivision)) == 85);
assert(double(J1.spectral_solve_count(k_subdivision)) == 255);
assert(double(J1.fem_refinement_used(k_subdivision)) == 1);
assert(double(J1.max_cr_refinement_level_used(k_subdivision)) >= 1);
assert(double(J1.max_cr_refinement_level_used(k_subdivision)) ...
    <= max_fem_refinement_level);
assert(double(J1.minimum_cr_mesh_size_used(k_subdivision)) ...
    == minimum_cr_mesh_size);
assert(double(J1.LG_fallback_node_count(k_subdivision)) ...
    <= double(J1.LG_CR_Liu_shift_node_count(k_subdivision)));
assert(contains(lower(string(J1.note(k_subdivision))),'subdivision'));
cell1_record = saved_records{find(cellfun( ...
    @(r) r.cell_id == 1,saved_records),1)};
assert(numel(cell1_record.leaves) == 64);
assert(all([cell1_record.leaves.depth] == 3));
assert(all(strcmp({cell1_record.leaves.method},'CR')));
assert(all([cell1_record.leaves.cr_refinement_level] == 2));
assert(all([cell1_record.leaves.lg_refinement_level] == -1));
assert(all([cell1_record.leaves.selected_attempt_ordinal] == 3));
assert(all([cell1_record.leaves.first_strict_attempt_ordinal] == 3));
assert(all(arrayfun(@(leaf) numel(leaf.attempts) == 3, ...
    cell1_record.leaves)));
fprintf(['cell 1 adaptive certificate: J1=%.17g, J2=%.17g, ', ...
    'leaves=%d, depth=%d, nodes=%d, solves=%d, hCR=%.17g\n'], ...
    cell1_record.J1_lower,cell1_record.J2_lower, ...
    cell1_record.leaf_count,cell1_record.max_leaf_depth, ...
    cell1_record.evaluated_node_count, ...
    cell1_record.spectral_solve_count, ...
    cell1_record.minimum_cr_mesh_size_used);

for cell_id = [138748,152387]
    k = find(double(J1.cell_id) == cell_id,1);
    assert(strcmp(char(J1.method(k)),'LG'));
    assert(isfinite(double(J1.LG_rho_lower(k))) ...
        && double(J1.LG_rho_lower(k)) > 0);
    assert(isfinite(double(J1.LG_ritz_upper(k))) ...
        && double(J1.LG_ritz_upper(k)) > 0);
    assert(double(J1.LG_rho_lower(k)) ...
        > double(J1.LG_ritz_upper(k)));
    assert(isfinite(double(J1.LG_separation_margin_lower(k))) ...
        && double(J1.LG_separation_margin_lower(k)) > 0);
    assert(isfinite(double(J1.LG_mu_upper(k))) ...
        && double(J1.LG_mu_upper(k)) < 0);
    assert(double(J1.LG_B_verified_spd(k)) == 1);
    shift_method = char(J1.LG_shift_method(k));
    assert(strcmp(shift_method,'CR-Liu'));
    assert(double(J1.LG_fallback_used(k)) == 0);
    fallback_reason = string(J1.LG_fallback_reason(k));
    assert(ismissing(fallback_reason) || strlength(fallback_reason) == 0);
end

k_cr = find(double(J1.cell_id) == 11731,1);
assert(strcmp(char(J1.method(k_cr)),'CR'));
assert(strcmp(char(J1.LG_shift_method(k_cr)),'not-applicable'));
assert(double(J1.LG_fallback_used(k_cr)) == 0);
assert(isnan(double(J1.LG_rho_lower(k_cr))));
assert(isnan(double(J1.LG_ritz_upper(k_cr))));
assert(isnan(double(J1.LG_separation_margin_lower(k_cr))));
assert(isnan(double(J1.LG_mu_upper(k_cr))));

json_text = fileread(json_path);
assert(~contains(json_text,project_root));
assert(~contains(json_text,output_dir));
assert(strcmp(manifest.hashes.J1_csv,ver10_file_sha256(j1_path)));
assert(strcmp(manifest.hashes.J2_csv,ver10_file_sha256(j2_path)));
assert(strcmp(manifest.hashes.geometric_coverage_json, ...
    ver10_file_sha256(coverage_json_path)));
assert(strcmp(manifest.hashes.mat,ver10_file_sha256(mat_path)));
assert(strcmp(manifest.hashes.config_mat, ...
    ver10_file_sha256(config_path)));
coverage_json = jsondecode(fileread(coverage_json_path));
assert(coverage_json.complete && coverage_json.rigorous);
assert(strcmp(coverage_json.hashes.proof_payload, ...
    manifest.geometric_coverage.hashes.proof_payload));

% All replay artifacts, including the MAT file that carries complete leaf
% and attempt histories, must be bound to the unified manifest.  Appending
% one byte to each artifact must make the publication validator fail before
% any copied generation can be accepted.
artifact_options = { ...
    'require_complete',false, ...
    'J1_file',j1_path,'J2_file',j2_path, ...
    'coverage_file',coverage_json_path, ...
    'mat_file',mat_path,'config_file',config_path, ...
    'unified_manifest_file',json_path};
binding = validate_omega_mid_unified_artifacts( ...
    manifest,project_root,artifact_options{:});
assert(binding.records_variable_count == numel(target_ids));
assert(strcmp(binding.leaf_attempt_evidence, ...
    'mat:records(:).leaves(:).attempts'));
tamper_dir = fullfile(output_dir,'artifact_tamper');
mkdir(tamper_dir);
tamper_sources = {j1_path,j2_path,coverage_json_path, ...
    mat_path,config_path};
tamper_options = {'J1_file','J2_file','coverage_file', ...
    'mat_file','config_file'};
for q = 1:numel(tamper_sources)
    [~,tamper_name,tamper_extension] = fileparts(tamper_sources{q});
    tampered = fullfile(tamper_dir, ...
        [tamper_name,'_tampered',tamper_extension]);
    [copied,message] = copyfile(tamper_sources{q},tampered,'f');
    assert(copied,message);
    local_append_tamper_byte(tampered);
    options = artifact_options;
    option_index = find(strcmp(options,tamper_options{q}),1);
    assert(~isempty(option_index));
    options{option_index+1} = tampered;
    local_assert_artifact_hash_rejected(@() ...
        validate_omega_mid_unified_artifacts( ...
            manifest,project_root,options{:}));
end

% A second invocation must use all four checkpoint records and perform no
% new cell solve.
resumed = run_omega_mid_unified_parallel( ...
    'mode','interval', ...
    'workers',0, ...
    'input_file',fullfile(project_root,'inputs','cell_def.csv'), ...
    'output_dir',output_dir, ...
    'run_name','targeted_regression', ...
    'chunk_size',2, ...
    'checkpoint_every',1, ...
    'max_subdivision_depth',max_subdivision_depth, ...
    'max_fem_refinement_level',max_fem_refinement_level, ...
    'minimum_cr_mesh_size',minimum_cr_mesh_size, ...
    'minimum_lg_mesh_size',minimum_lg_mesh_size, ...
    'cell_ids',target_ids, ...
    'resume',true, ...
    'retry_failed',true);
assert(resumed.summary.num_resumed_records == 4);
assert(resumed.summary.num_computed_this_invocation == 0);
assert(resumed.summary.exact_selected_id_coverage);

% A boolean-only checkpoint must not be trusted.  Corrupt a persisted leaf
% certificate for former-Inf cell 138748 while leaving the parent minimum
% and all flags true.  Strict reload must validate every leaf and recompute
% exactly that parent cell even when retry_failed=false: that option may
% retain a genuinely failed record, but may never trust a claimed success.
checkpoint_files = dir(fullfile( ...
    run_dir,'checkpoints','chunk_*.mat'));
checkpoint_path = '';
payload = struct();
k_bad = [];
for q = 1:numel(checkpoint_files)
    candidate_path = fullfile( ...
        checkpoint_files(q).folder,checkpoint_files(q).name);
    candidate = load(candidate_path);
    stored_ids = cellfun(@(r) r.cell_id,candidate.records);
    candidate_bad = find(stored_ids == 138748,1);
    if ~isempty(candidate_bad)
        checkpoint_path = candidate_path;
        payload = candidate;
        k_bad = candidate_bad;
        break;
    end
end
assert(~isempty(checkpoint_path));
assert(~isempty(k_bad));
assert(~isempty(payload.records{k_bad}.leaves));
selected_attempt = ...
    payload.records{k_bad}.leaves(1).selected_attempt_ordinal;
assert(selected_attempt >= 1);
payload.records{k_bad}.leaves(1).attempts( ...
    selected_attempt).J1_lower = NaN;
save(checkpoint_path,'-struct','payload','-v7');

repaired = run_omega_mid_unified_parallel( ...
    'mode','interval', ...
    'workers',0, ...
    'input_file',fullfile(project_root,'inputs','cell_def.csv'), ...
    'output_dir',output_dir, ...
    'run_name','targeted_regression', ...
    'chunk_size',2, ...
    'checkpoint_every',1, ...
    'max_subdivision_depth',max_subdivision_depth, ...
    'max_fem_refinement_level',max_fem_refinement_level, ...
    'minimum_cr_mesh_size',minimum_cr_mesh_size, ...
    'minimum_lg_mesh_size',minimum_lg_mesh_size, ...
    'cell_ids',target_ids, ...
    'resume',true, ...
    'retry_failed',false);
assert(repaired.summary.num_resumed_records == 3);
assert(repaired.summary.num_computed_this_invocation == 1);
assert(repaired.summary.num_strict_records == 4);
J1_repaired = readtable(j1_path,'TextType','string');
assert(all(isfinite(double(J1_repaired.J_lower)) ...
    & double(J1_repaired.J_lower) > 0));

% Terminal leaves are fail-closed.  With no geometric or FEM refinement,
% cell 1 cannot close; the driver must preserve the failed CR/LG attempts,
% mark the parent invalid, and publish verified=0 with no numerical bound.
failed = run_omega_mid_unified_parallel( ...
    'mode','interval', ...
    'workers',0, ...
    'input_file',fullfile(project_root,'inputs','cell_def.csv'), ...
    'output_dir',output_dir, ...
    'run_name','terminal_fail_closed', ...
    'chunk_size',1, ...
    'checkpoint_every',1, ...
    'max_subdivision_depth',0, ...
    'max_fem_refinement_level',0, ...
    'minimum_cr_mesh_size',minimum_cr_mesh_size, ...
    'minimum_lg_mesh_size',minimum_lg_mesh_size, ...
    'cell_ids',1, ...
    'resume',true, ...
    'retry_failed',true);
assert(failed.summary.num_strict_records == 0);
assert(failed.summary.num_invalid_results == 1);
assert(failed.summary.num_errors == 0);
assert(~failed.summary.exact_full_input_id_coverage);
assert(~failed.complete_certificate);
failed_dir = fullfile(output_dir,'terminal_fail_closed');
failed_payload = load(fullfile( ...
    failed_dir,'omega_mid_unified_results.mat'),'records');
failed_record = failed_payload.records{1};
assert(strcmp(failed_record.status,'invalid_result'));
assert(~failed_record.strict_checks_passed);
assert(isnan(failed_record.J1_lower));
assert(isnan(failed_record.J2_lower));
assert(failed_record.leaf_count == 1);
assert(numel(failed_record.leaves(1).attempts) == 2);
assert(failed_record.leaves(1).first_strict_attempt_ordinal == 0);
assert(~failed_record.leaves(1).strict_checks_passed);
failed_csv = readtable(fullfile( ...
    failed_dir,'J1_OmegaMid.csv'),'TextType','string');
assert(double(failed_csv.verified(1)) == 0);
assert(strcmp(char(failed_csv.status(1)),'invalid_result'));
assert(isnan(double(failed_csv.J_lower(1))));
fprintf(['terminal fail-closed certificate: status=%s, ', ...
    'strict=%d, attempts=%d, J1=%.17g, J2=%.17g\n'], ...
    failed_record.status,failed_record.strict_checks_passed, ...
    failed_record.leaves(1).eigenvalue_attempts, ...
    failed_record.J1_lower,failed_record.J2_lower);

fprintf(['test_omega_mid_unified_driver: PASS ', ...
    '(CR=1, LG=2, subdivision=1, resumed=4, ', ...
    'corrupt-leaf checkpoint repaired=1, ', ...
    'artifact tamper cases=5, terminal fail-closed=1)\n']);
clear cleanup
local_remove_tree(output_dir);
end


function local_append_tamper_byte(filename)
fid = fopen(filename,'ab');
if fid < 0
    error('test_omega_mid_unified_driver:TamperOpenFailed', ...
        'Could not open tamper fixture %s.',filename);
end
cleanup = onCleanup(@() fclose(fid));
fwrite(fid,uint8(10),'uint8');
clear cleanup
end


function local_assert_artifact_hash_rejected(callback)
caught = false;
try
    callback();
catch ME
    if ~strcmp(ME.identifier, ...
            ['validate_omega_mid_unified_artifacts:', ...
             'ArtifactHashMismatch'])
        rethrow(ME);
    end
    caught = true;
end
assert(caught);
end


function local_remove_tree(pathname)
if exist(pathname,'dir')
    rmdir(pathname,'s');
end
end
