function test_runtime_metadata_stability()
%TEST_RUNTIME_METADATA_STABILITY Repeated INTLAB saves preserve provenance.

if ~strcmp(getenv('RUN_INTLAB_SMOKE'),'1')
    fprintf('Runtime metadata stability test skipped without INTLAB.\n');
    return
end
project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(project_root,'scripts_run'));
snapshots = cell(3,1);
for k = 1:3
    omega_up_all_prepare_worker(project_root,'interval');
    snapshots{k} = ver10_runtime_metadata();
    assert(strcmp(snapshots{k}.intlab_tree_hash_policy,[ ...
        'recursive-path-size-content-sha256-v2-with-exact-intlab-data-', ...
        'mat-v5-created-on-timestamp-normalization']));
    records = snapshots{k}.intlab_data_mat_records;
    assert(numel(records) == 1);
    assert(strcmp(records.normalization, ...
        'mat-v5-created-on-24-character-timestamp-zeroed-v1'));
    assert(records.bytes >= 128);
    assert(~isempty(regexp(records.normalized_sha256, ...
        '^[0-9a-f]{64}$','once')));
end
assert(strcmp(jsonencode(snapshots{1}),jsonencode(snapshots{2})) ...
    && strcmp(jsonencode(snapshots{1}),jsonencode(snapshots{3})));
fprintf(['test_runtime_metadata_stability: PASS ', ...
    '(three repeated INTLAB saves, identical normalized runtime)\n']);
end
