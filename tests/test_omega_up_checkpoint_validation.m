function test_omega_up_checkpoint_validation()
%TEST_OMEGA_UP_CHECKPOINT_VALIDATION Reject mismatched/corrupt resume data.
%
% The selected rectangles are rigorously outside the unit disk, so they
% exercise the full checkpoint/merge/manifest path without an FEM solve.
% After the first run, one checkpoint receives the wrong task id and the
% other receives an impossible +Inf functional bound.  Both must be
% rejected and recomputed; a subsequent clean resume must load both.

project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(project_root,'scripts_run'));

output_dir = tempname;
mkdir(output_dir);
cleanup = onCleanup(@() local_remove_tree(output_dir));

task_ids = [29525,29647]; % ix=243,244 and iy=1 for Nx=244, Ny=122
common = { ...
    'mode','double', ...
    'workers',0, ...
    'output_dir',output_dir, ...
    'run_name','checkpoint_validation', ...
    'eps_up','0.122', ...
    'Nx',244, ...
    'Ny',122, ...
    'Ny_axis',488, ...
    'functional_scope','split', ...
    'jup_y_min','0.85', ...
    'task_ids',task_ids, ...
    'resume',true, ...
    'retry_failed',true};

first = run_omega_up_all_residual_parallel(common{:});
assert(first.summary.num_resumed == 0);
assert(first.summary.exact_selected_task_id_coverage);
assert(~first.summary.exact_full_task_id_coverage);
assert(~first.summary.complete_certificate);
assert(first.runtime_rechecked_after_run);
assert(ischar(first.runtime_sha256) ...
    && ~isempty(regexp(first.runtime_sha256, ...
        '^[0-9a-f]{64}$','once')));
assert(strcmp(first.config.runtime_sha256,first.runtime_sha256));

checkpoint_dir = fullfile( ...
    output_dir,'checkpoint_validation','checkpoints');
id_checkpoint = fullfile(checkpoint_dir,sprintf( ...
    'task_%06d.mat',task_ids(1)));
inf_checkpoint = fullfile(checkpoint_dir,sprintf( ...
    'task_%06d.mat',task_ids(2)));

payload = load(id_checkpoint,'record','fingerprint');
record = payload.record;
fingerprint = payload.fingerprint;
record.task_id = task_ids(2);
save(id_checkpoint,'record','fingerprint','-v7');

payload = load(inf_checkpoint,'record','fingerprint');
record = payload.record;
fingerprint = payload.fingerprint;
record.J1_lower = Inf;
save(inf_checkpoint,'record','fingerprint','-v7');

repaired = run_omega_up_all_residual_parallel(common{:});
assert(repaired.summary.num_resumed == 0);
assert(repaired.summary.exact_selected_task_id_coverage);
assert(repaired.summary.num_duplicate_task_ids == 0);
assert(repaired.summary.num_missing_task_ids == 0);
assert(repaired.summary.num_unexpected_task_ids == 0);

payload = load(id_checkpoint,'record');
assert(payload.record.task_id == task_ids(1));
payload = load(inf_checkpoint,'record');
assert(payload.record.task_id == task_ids(2));
assert(isnan(payload.record.J1_lower));

resumed = run_omega_up_all_residual_parallel(common{:});
assert(resumed.summary.num_resumed == numel(task_ids));
assert(resumed.summary.exact_selected_task_id_coverage);
assert(resumed.runtime_rechecked_after_run);
assert(strcmp(resumed.config.runtime_sha256,resumed.runtime_sha256));

fprintf(['test_omega_up_checkpoint_validation: PASS ', ...
    '(bad id and +Inf checkpoints rejected)\n']);
clear cleanup
local_remove_tree(output_dir);
end


function local_remove_tree(pathname)
if exist(pathname,'dir')
    rmdir(pathname,'s');
end
end
