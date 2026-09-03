%RUN_OMEGA_UP_SPECTRAL_ATLAS  Rigorous reusable CR--LG endpoint atlas.
%
% Required:
%   VER10_SPECTRAL_ATLAS=/absolute/output/omega_up_spectral_atlas.mat
%
% Optional:
%   VER10_ATLAS_NX       number of x anchors (default 9)
%   VER10_ATLAS_NY       number of y anchors (default 5)
%   VER10_ATLAS_WORKERS  process workers (default 20)

project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(project_root,'scripts_run'));
output_path = getenv('VER10_SPECTRAL_ATLAS');
if isempty(output_path)
    error('run_omega_up_spectral_atlas:MissingOutput', ...
        'Set VER10_SPECTRAL_ATLAS to an absolute MAT-file path.');
end
if exist(output_path,'file')
    error('run_omega_up_spectral_atlas:OutputExists', ...
        'The requested atlas file already exists: %s',output_path);
end

Nx_anchor = local_environment_integer('VER10_ATLAS_NX',9);
Ny_anchor = local_environment_integer('VER10_ATLAS_NY',5);
workers = local_environment_integer('VER10_ATLAS_WORKERS',20);
if Nx_anchor < 1 || Ny_anchor < 1
    error('run_omega_up_spectral_atlas:BadAnchorCount', ...
        'VER10_ATLAS_NX and VER10_ATLAS_NY must be positive.');
end
N_LG = 16;
N_rho = 64;
bound_order = 2;

source_before = ver10_source_provenance(project_root);
if ~source_before.git_status_available || source_before.git_dirty ...
        || isempty(source_before.git_commit)
    error('run_omega_up_spectral_atlas:DirtySource', ...
        'The rigorous atlas must start from a clean committed source tree.');
end
omega_up_all_prepare_worker(project_root,'interval');
global INTERVAL_MODE
if ~logical(INTERVAL_MODE)
    error('run_omega_up_spectral_atlas:IntervalInitFailed', ...
        'INTLAB interval mode was not enabled.');
end

% Freeze the complete client runtime before forming any interval endpoint.
% Every worker is compared with this exact JSON record after initialization
% and again after its spectral work has completed.
runtime = ver10_runtime_metadata();
runtime_payload = jsonencode(runtime);
runtime_sha256 = local_sha256_bytes(uint8(runtime_payload));

eps_up = I_intval('0.122');
x_left = I_intval('0.5');
x_right = x_left+2*eps_up;
y_top = sqrt(I_intval('3'))/2;
y_bottom = y_top-eps_up;
x_anchors = local_anchor_vector(x_left,x_right,Nx_anchor);
y_anchors = local_anchor_vector(y_bottom,y_top,Ny_anchor);
num_anchors = Nx_anchor*Ny_anchor;
certificates_linear = cell(num_anchors,1);

started = tic;
worker_runtime_before = local_serial_worker_report( ...
    runtime_sha256,'spectral-atlas-before-dispatch');
worker_runtime_after = local_serial_worker_report( ...
    runtime_sha256,'spectral-atlas-after-dispatch');
if workers > 0
    pool = ver10_ensure_local_pool(workers);
    future = parfevalOnAll( ...
        pool,@omega_up_all_prepare_worker,0,project_root,'interval');
    wait(future);
    fetchOutputs(future);
    worker_runtime_before = ver10_assert_worker_runtime( ...
        pool,runtime,runtime_sha256,'spectral-atlas-before-dispatch');
    parfor k = 1:num_anchors
        [iy,ix] = ind2sub([Ny_anchor,Nx_anchor],k);
        triangle = I_intval( ...
            [0,0,1,0,x_anchors(ix),y_anchors(iy)]);
        certificates_linear{k} = ...
            triangle_spectral_certificate_cr_lg( ...
                triangle,N_LG,N_rho,bound_order);
    end
    worker_runtime_after = ver10_assert_worker_runtime( ...
        pool,runtime,runtime_sha256,'spectral-atlas-after-dispatch');
else
    for k = 1:num_anchors
        [iy,ix] = ind2sub([Ny_anchor,Nx_anchor],k);
        triangle = I_intval( ...
            [0,0,1,0,x_anchors(ix),y_anchors(iy)]);
        certificates_linear{k} = ...
            triangle_spectral_certificate_cr_lg( ...
                triangle,N_LG,N_rho,bound_order);
        fprintf('spectral atlas: %d/%d anchors complete\n', ...
            k,num_anchors);
    end
end
wall_seconds = toc(started);
certificate_rigorous = cellfun( ...
    @(c) isfield(c,'rigorous') ...
        && islogical(c.rigorous) ...
        && isscalar(c.rigorous) ...
        && c.rigorous, ...
    certificates_linear);
if ~all(certificate_rigorous)
    error('run_omega_up_spectral_atlas:UnverifiedCertificate', ...
        'At least one anchor did not return a rigorous certificate.');
end

ver10_assert_runtime_snapshot( ...
    runtime,runtime_sha256,'spectral-atlas-client-after-run');
source_after = ver10_assert_source_unchanged(project_root,source_before);
atlas = struct();
atlas.schema = ...
    'lowerboundsineq.omega-up-spectral-cr-lg-atlas.v1';
atlas.method = ...
    'nearest-anchor-affine-transport-of-LG-Ritz-CR-endpoints';
atlas.rigorous = logical(all(certificate_rigorous));
atlas.eps_up = eps_up;
atlas.x_domain = I_infsup(I_inf(x_left),I_sup(x_right));
atlas.y_domain = I_infsup(I_inf(y_bottom),I_sup(y_top));
atlas.x_domain_lower = I_inf(x_left);
atlas.x_domain_upper = I_sup(x_right);
atlas.y_domain_lower = I_inf(y_bottom);
atlas.y_domain_upper = I_sup(y_top);
atlas.x_anchors = x_anchors;
atlas.y_anchors = y_anchors;
atlas.certificates = ...
    reshape(certificates_linear,[Ny_anchor,Nx_anchor]);
atlas.Nx_anchor = Nx_anchor;
atlas.Ny_anchor = Ny_anchor;
atlas.N_LG = N_LG;
atlas.N_rho = N_rho;
atlas.bound_order = bound_order;
atlas.num_certificates = num_anchors;
atlas.wall_seconds = wall_seconds;
atlas.source = source_after;
atlas.runtime = runtime;
atlas.runtime_sha256 = runtime_sha256;
atlas.runtime_rechecked_after_run = true;
atlas.parallel_workers = workers;
atlas.worker_runtime_validation = struct( ...
    'policy','full-runtime-json-all-workers-before-and-after-dispatch-v1', ...
    'before_dispatch',worker_runtime_before, ...
    'after_dispatch',worker_runtime_after, ...
    'all_computation_workers_checked',true, ...
    'runtime_sha256',runtime_sha256);

output_directory = fileparts(output_path);
if ~isempty(output_directory) && ~exist(output_directory,'dir')
    mkdir(output_directory);
end
temporary_path = sprintf('%s.tmp.%d',output_path,feature('getpid'));
cleanup = onCleanup(@() local_remove_temporary(temporary_path));
save(temporary_path,'atlas','-v7.3');
[ok,message] = movefile(temporary_path,output_path);
if ~ok
    error('run_omega_up_spectral_atlas:MoveFailed','%s',message);
end
clear cleanup

fprintf('rigorous CR--LG spectral atlas complete: %s\n',output_path);
fprintf('anchors: %d, wall seconds: %.3f\n',num_anchors,wall_seconds);
fprintf('SHA-256: %s\n',ver10_file_sha256(output_path));


function value = local_environment_integer(name,default_value)
text = getenv(name);
if isempty(text)
    value = default_value;
else
    value = str2double(text);
end
if ~isscalar(value) || ~isfinite(value) ...
        || value < 0 || value ~= floor(value)
    error('run_omega_up_spectral_atlas:BadEnvironment', ...
        '%s must be a nonnegative integer.',name);
end
end


function anchors = local_anchor_vector(left,right,count)
if count == 1
    anchors = I_mid((left+right)/2);
else
    anchors = linspace(I_mid(left),I_mid(right),count);
end
end


function local_remove_temporary(path)
if exist(path,'file')
    delete(path);
end
end


function report = local_serial_worker_report(runtime_sha256,phase)
report = struct( ...
    'schema','lowerboundsineq.worker-runtime-validation.v1', ...
    'phase',phase, ...
    'runtime_sha256',runtime_sha256, ...
    'worker_count',0, ...
    'all_workers_returned',true, ...
    'all_workers_full_runtime_json_match',true);
end


function digest = local_sha256_bytes(bytes)
try
    engine = javaMethod( ...
        'getInstance','java.security.MessageDigest','SHA-256');
catch ME
    error('run_omega_up_spectral_atlas:NoSHA256', ...
        'Cannot initialize SHA-256: %s',ME.message);
end
bytes = uint8(bytes(:));
if ~isempty(bytes)
    engine.update(typecast(bytes,'int8'));
end
raw = typecast(engine.digest(),'uint8');
digest = lower(reshape(dec2hex(raw,2).',1,[]));
end
