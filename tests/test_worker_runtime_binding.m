function test_worker_runtime_binding()
%TEST_WORKER_RUNTIME_BINDING Exercise all-worker full-runtime validation.

project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(project_root,'scripts_run'));
omega_up_all_prepare_worker(project_root,'double');

trusted = ver10_runtime_metadata();
trusted_sha256 = local_sha256_bytes(uint8(jsonencode(trusted)));
pool = ver10_ensure_local_pool(2);
cleanup = onCleanup(@() local_delete_pool(pool));
initialization = parfevalOnAll( ...
    pool,@omega_up_all_prepare_worker,0,project_root,'double');
wait(initialization);
fetchOutputs(initialization);

before = ver10_assert_worker_runtime( ...
    pool,trusted,trusted_sha256,'worker-runtime-test-before');
after = ver10_assert_worker_runtime( ...
    pool,trusted,trusted_sha256,'worker-runtime-test-after');
assert(before.worker_count == 2 && after.worker_count == 2);
assert(strcmp(before.broadcast_output_class,'struct') ...
    && before.broadcast_output_count == 2, ...
    'Unexpected R2023b FevalOnAllFuture broadcast output shape.');
assert(before.all_workers_returned && after.all_workers_returned);
assert(before.all_workers_full_runtime_json_match ...
    && after.all_workers_full_runtime_json_match);
record = struct( ...
    'policy','full-runtime-json-all-workers-before-and-after-dispatch-v1', ...
    'before_dispatch',before, ...
    'after_dispatch',after, ...
    'all_computation_workers_checked',true, ...
    'runtime_sha256',trusted_sha256);
ver10_validate_worker_runtime_record( ...
    record,trusted_sha256,2,'synthetic worker record');

bad_record = record;
bad_record.all_computation_workers_checked = 2;
rejected_bad_boolean = false;
try
    ver10_validate_worker_runtime_record( ...
        bad_record,trusted_sha256,2,'tampered worker record');
catch ME
    rejected_bad_boolean = strcmp(ME.identifier, ...
        'ver10_validate_worker_runtime_record:InvalidRecord');
end
assert(rejected_bad_boolean, ...
    'A numeric certificate-critical worker flag was accepted.');

tampered = trusted;
tampered.gmsh_version = [char(tampered.gmsh_version),'-tampered'];
tampered_sha256 = local_sha256_bytes(uint8(jsonencode(tampered)));
rejected = false;
try
    ver10_assert_worker_runtime( ...
        pool,tampered,tampered_sha256,'worker-runtime-test-tamper');
catch ME
    fprintf('tamper rejection wrapper: %s | %s\n', ...
        ME.identifier,ME.message);
    rejected = local_has_runtime_mismatch(ME);
end
assert(rejected,'A worker accepted a tampered trusted runtime.');

fprintf([ ...
    'test_worker_runtime_binding: PASS (2 workers, broadcast=%s[%d], ', ...
    'tamper rejected)\n'],before.broadcast_output_class, ...
    before.broadcast_output_count);
clear cleanup
local_delete_pool(pool);
end


function local_delete_pool(pool)
if ~isempty(pool) && isvalid(pool)
    delete(pool);
end
end


function digest = local_sha256_bytes(bytes)
engine = javaMethod( ...
    'getInstance','java.security.MessageDigest','SHA-256');
bytes = uint8(bytes(:));
if ~isempty(bytes)
    engine.update(typecast(bytes,'int8'));
end
raw = typecast(engine.digest(),'uint8');
digest = lower(reshape(dec2hex(raw,2).',1,[]));
end


function tf = local_has_runtime_mismatch(exception)
tf = strcmp(exception.identifier, ...
        'ver10_assert_runtime_snapshot:RuntimeMismatch') ...
    || contains(exception.message, ...
        'process runtime differs from the frozen client runtime');
if tf
    return;
end
causes = exception.cause;
for k = 1:numel(causes)
    if local_has_runtime_mismatch(causes{k})
        tf = true;
        return;
    end
end
end
