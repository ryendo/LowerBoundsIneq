function report = ver10_assert_runtime_snapshot( ...
    trusted_runtime,expected_sha256,phase)
%VER10_ASSERT_RUNTIME_SNAPSHOT Match this MATLAB process to a frozen runtime.

if nargin < 3
    phase = 'unspecified';
end
phase = char(phase);
if ~isstruct(trusted_runtime) || ~isscalar(trusted_runtime)
    error('ver10_assert_runtime_snapshot:BadTrustedRuntime', ...
        'The trusted runtime must be a scalar structure.');
end
expected_sha256 = char(expected_sha256);
if isempty(regexp(expected_sha256,'^[0-9a-f]{64}$','once'))
    error('ver10_assert_runtime_snapshot:BadExpectedHash', ...
        'The expected runtime SHA-256 is malformed.');
end

trusted_payload = jsonencode(trusted_runtime);
trusted_sha256 = local_sha256_bytes(uint8(trusted_payload));
if ~strcmp(trusted_sha256,expected_sha256)
    error('ver10_assert_runtime_snapshot:TrustedHashMismatch', ...
        'The trusted runtime structure does not match its SHA-256.');
end

current_runtime = ver10_runtime_metadata();
current_payload = jsonencode(current_runtime);
current_sha256 = local_sha256_bytes(uint8(current_payload));
if ~strcmp(current_payload,trusted_payload) ...
        || ~strcmp(current_sha256,expected_sha256)
    error('ver10_assert_runtime_snapshot:RuntimeMismatch', ...
        ['The %s process runtime differs from the frozen client runtime. ', ...
         'The interval proof cannot mix MATLAB, INTLAB, Gmsh, operating-', ...
         'system, or hardware states.'],phase);
end

report = struct( ...
    'schema','lowerboundsineq.runtime-snapshot-validation.v1', ...
    'phase',phase, ...
    'runtime_sha256',current_sha256, ...
    'full_runtime_json_match',true);
end


function digest = local_sha256_bytes(bytes)
try
    engine = javaMethod( ...
        'getInstance','java.security.MessageDigest','SHA-256');
catch ME
    error('ver10_assert_runtime_snapshot:NoSHA256', ...
        'Cannot initialize SHA-256: %s',ME.message);
end
bytes = uint8(bytes(:));
if ~isempty(bytes)
    engine.update(typecast(bytes,'int8'));
end
raw = typecast(engine.digest(),'uint8');
digest = lower(reshape(dec2hex(raw,2).',1,[]));
end
