function manifest = run_omega_up_all_residual_parallel(varargin)
%RUN_OMEGA_UP_ALL_RESIDUAL_PARALLEL  Resumable 20-worker Omega_up proof.
%
% The default rigorous run uses dx=dy=0.00025 rectangles and a finer
% dy_axis=0.0000625 symmetry-axis grid (Ny_axis=1952), twenty workers, and the
% degree-11 global-Bernstein strong-residual Hessian estimator.  It uses no
% mesh or H(div) flux reconstruction in the Hessian step.  A reusable
% coarse CR--Liu/Lehmann--Goerisch atlas supplies the three spectral
% endpoints.  Each rectangle invokes the
% x-direction estimator once and each symmetry-axis interval invokes the
% y-direction estimator once.  In the default 'split' scope the same result
% tests J1 and J2 on all of Omega_up and Siudeja's upper functional JUP only
% on the local band y>=0.85.  Out-of-scope functionals are recorded
% explicitly and never counted as failures.
%
% Examples:
%   % Full liulabhpc certificate (defaults shown explicitly):
%   run_omega_up_all_residual_parallel( ...
%       'mode','interval','workers',20,'dx','0.00025','dy','0.00025', ...
%       'Ny_axis',1952,'functional_scope','split','jup_y_min','0.85');
%
%   % Small exploratory serial run:
%   run_omega_up_all_residual_parallel( ...
%       'mode','double','workers',0,'max_tasks',3, ...
%       'run_name','local_smoke');
%
% Checkpointing:
%   Every task owns one MAT checkpoint below RUN_DIR/checkpoints.  A killed
%   job can therefore be resumed without trusting a partially written
%   shared CSV.  On every completed invocation the checkpoints are merged
%   deterministically into CSV, MAT, and a JSON manifest.

project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(project_root,'scripts_run'));
source_provenance = ver10_source_provenance(project_root);

p = inputParser;
p.CaseSensitive = true;
if isprop(p,'PartialMatching')
    p.PartialMatching = false;
end
addParameter(p,'mode','interval',@(x) ischar(x) || isstring(x));
addParameter(p,'eps_up','0.122',@local_positive_scalar_input);
addParameter(p,'dx','0.00025',@local_positive_scalar_input);
addParameter(p,'dy','0.00025',@local_positive_scalar_input);
addParameter(p,'functional_scope','split',@local_scope_name_input);
addParameter(p,'jup_y_min','0.85',@local_positive_scalar_input);
addParameter(p,'Nx',[],@local_optional_positive_integer);
addParameter(p,'Ny',[],@local_optional_positive_integer);
addParameter(p,'Ny_axis',1952,@local_optional_positive_integer);
addParameter(p,'workers',20,@local_nonnegative_integer);
addParameter(p,'mesh_params',struct(),@isstruct);
addParameter(p,'output_dir', ...
    fullfile(project_root,'results','omega_up_all_residual'), ...
    @(x) ischar(x) || isstring(x));
addParameter(p,'run_name','',@(x) ischar(x) || isstring(x));
addParameter(p,'resume',true,@(x) islogical(x) || isnumeric(x));
addParameter(p,'retry_failed',false,@(x) islogical(x) || isnumeric(x));
addParameter(p,'max_tasks',Inf, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p,'task_ids',[],@local_optional_task_ids);
parse(p,varargin{:});
opt = p.Results;

mode = lower(char(opt.mode));
if ~any(strcmp(mode,{'interval','double'}))
    error('run_omega_up_all_residual_parallel:BadMode', ...
        'mode must be ''interval'' or ''double''.');
end

% Initialize INTLAB before loading a production atlas: its MAT payload
% contains nested intval objects and must never be deserialized before the
% class is available.
omega_up_all_prepare_worker(project_root,mode);
global INTERVAL_MODE
rigorous = logical(INTERVAL_MODE);
if strcmp(mode,'interval') && ~rigorous
    error('run_omega_up_all_residual_parallel:IntervalInitFailed', ...
        'The requested interval runtime did not enable INTERVAL_MODE.');
end
runtime_metadata = ver10_runtime_metadata();
runtime_payload = jsonencode(runtime_metadata);
runtime_sha256 = local_sha256_bytes(uint8(runtime_payload));

output_dir = char(opt.output_dir);
resume = logical(opt.resume);
retry_failed = logical(opt.retry_failed);
mesh_params = local_default_mesh_parameters(opt.mesh_params);
if isfield(mesh_params,'spectral_atlas') ...
        && ~isempty(mesh_params.spectral_atlas)
    error('run_omega_up_all_residual_parallel:InMemoryAtlasUnsupported', ...
        ['Resumable runs require spectral_atlas_path; an in-memory ', ...
         'atlas cannot be fingerprinted safely.']);
end
spectral_atlas_hash = '';
spectral_atlas_metadata = struct();
if isfield(mesh_params,'spectral_atlas_path') ...
        && ~isempty(mesh_params.spectral_atlas_path)
    spectral_atlas_path = char(mesh_params.spectral_atlas_path);
    if ~exist(spectral_atlas_path,'file')
        error('run_omega_up_all_residual_parallel:MissingSpectralAtlas', ...
            'The spectral-atlas MAT file does not exist.');
    end
    spectral_atlas_hash = ver10_file_sha256(spectral_atlas_path);
    spectral_atlas_metadata = local_validate_spectral_atlas( ...
        spectral_atlas_path,spectral_atlas_hash, ...
        source_provenance.git_commit,strcmp(mode,'interval'));
end
mesh_params.spectral_atlas_sha256 = spectral_atlas_hash;

epsI = local_decimal_value(opt.eps_up);
dxI = local_decimal_value(opt.dx);
dyI = local_decimal_value(opt.dy);
jup_y_minI = local_decimal_value(opt.jup_y_min);
if ~(I_inf(epsI) > 0) || ~(I_inf(dxI) > 0) || ~(I_inf(dyI) > 0)
    error('run_omega_up_all_residual_parallel:BadGrid', ...
        'eps_up, dx, and dy must be strictly positive.');
end
y_top = sqrt(I_intval('3'))/2;
if I_inf(jup_y_minI) > I_sup(y_top)
    error('run_omega_up_all_residual_parallel:EmptyJupBand', ...
        'jup_y_min lies above the top of Omega_up.');
end
scope = local_scope_config( ...
    lower(char(opt.functional_scope)),jup_y_minI);

if isempty(opt.Nx)
    Nx = ceil(I_sup(2*epsI/dxI));
else
    Nx = double(opt.Nx);
end
if isempty(opt.Ny)
    Ny = ceil(I_sup(epsI/dyI));
else
    Ny = double(opt.Ny);
end
if isempty(opt.Ny_axis)
    Ny_axis = ceil(I_sup(epsI/dyI));
else
    Ny_axis = double(opt.Ny_axis);
end

[cells,grid] = local_build_tasks(epsI,Nx,Ny,Ny_axis,scope);
domain_coverage = local_domain_coverage(cells,grid,scope);
full_task_count = numel(cells);
if ~isempty(opt.task_ids) && isfinite(opt.max_tasks)
    error('run_omega_up_all_residual_parallel:ConflictingTaskSelection', ...
        'Use either task_ids or max_tasks, not both.');
end
if ~isempty(opt.task_ids)
    selected_ids = unique(double(opt.task_ids(:).'),'stable');
    if any(selected_ids > full_task_count)
        error('run_omega_up_all_residual_parallel:BadTaskIds', ...
            'Every task id must be at most %d.',full_task_count);
    end
    cells = cells(selected_ids);
    selection_mode = 'explicit-task-ids';
    selection_signature = ['ids-',sprintf('%d_',selected_ids)];
elseif isfinite(opt.max_tasks)
    cells = cells(1:min(full_task_count,floor(opt.max_tasks)));
    selection_mode = 'prefix';
    selection_signature = sprintf('prefix-1-%d',numel(cells));
else
    selection_mode = 'all';
    selection_signature = 'all';
end
coverage_complete = numel(cells) == full_task_count ...
    && isequal([cells.id],1:full_task_count);
if strcmp(mode,'interval') && coverage_complete ...
        && strcmpi(char(mesh_params.estimator_backend),'bernstein-strong') ...
        && isempty(spectral_atlas_hash)
    error('run_omega_up_all_residual_parallel:SpectralAtlasRequired', ...
        ['A complete Bernstein interval run requires ', ...
         'mesh_params.spectral_atlas_path.']);
end
if strcmp(mode,'interval') && coverage_complete ...
        && (~source_provenance.git_status_available ...
            || source_provenance.git_dirty ...
            || isempty(source_provenance.git_commit))
    error('run_omega_up_all_residual_parallel:DirtyProductionTree', ...
        ['A complete interval run must start from a clean Git tree. ', ...
         'Commit the implementation before creating production checkpoints.']);
end

config = struct();
config.schema = 'lowerboundsineq.omega_up_all_residual.config.v3';
config.estimator_implementation = ...
    'signed-strong-residual-Bernstein-CRLG-atlas-split-v7';
config.mode = mode;
config.rigorous = rigorous;
config.eps_up = I_mid(epsI);
config.eps_up_lower = I_inf(epsI);
config.eps_up_upper = I_sup(epsI);
config.dx_requested = I_mid(dxI);
config.dx_requested_lower = I_inf(dxI);
config.dx_requested_upper = I_sup(dxI);
config.dy_requested = I_mid(dyI);
config.dy_requested_lower = I_inf(dyI);
config.dy_requested_upper = I_sup(dyI);
config.Nx = Nx;
config.Ny = Ny;
config.Ny_axis = Ny_axis;
config.dx_actual = 2*I_mid(epsI)/Nx;
config.dy_actual = I_mid(epsI)/Ny;
config.dy_axis_actual = I_mid(epsI)/Ny_axis;
config.functional_scope = scope;
config.mesh_params = mesh_params;
config.spectral_atlas_sha256 = spectral_atlas_hash;
if isempty(fieldnames(spectral_atlas_metadata))
    config.spectral_atlas_source_commit = '';
    config.spectral_atlas_schema = '';
else
    config.spectral_atlas_source_commit = ...
        spectral_atlas_metadata.source.git_commit;
    config.spectral_atlas_schema = spectral_atlas_metadata.schema;
end
config.workers = double(opt.workers);
config.full_task_count = full_task_count;
config.selected_task_count = numel(cells);
config.selection_mode = selection_mode;
config.selection_signature = selection_signature;
if strcmp(selection_mode,'explicit-task-ids')
    config.selected_task_ids = [cells.id];
else
    config.selected_task_ids = [];
end
config.coverage_complete = coverage_complete;
config.domain_coverage = domain_coverage;
config.runtime = runtime_metadata;
config.runtime_sha256 = runtime_sha256;
config.worker_runtime_validation_policy = ...
    'full-runtime-json-all-workers-before-and-after-dispatch-v1';
config.source_git_commit = source_provenance.git_commit;
config.source_git_dirty = source_provenance.git_dirty;
config.domain = ['Omega_up subset of {1/2<=x<=1/2+2 eps, ', ...
    'sqrt(3)/2-eps<=y<=sqrt(3)/2, x^2+y^2<=1}'];
fingerprint = local_config_fingerprint(config);

if isempty(char(opt.run_name))
    run_name = local_default_run_name(config);
else
    run_name = local_safe_name(char(opt.run_name));
end
run_dir = fullfile(output_dir,run_name);
checkpoint_dir = fullfile(run_dir,'checkpoints');
if ~exist(checkpoint_dir,'dir')
    mkdir(checkpoint_dir);
elseif ~resume
    existing = dir(fullfile(checkpoint_dir,'task_*.mat'));
    if ~isempty(existing)
        error('run_omega_up_all_residual_parallel:ExistingCheckpoints', ...
            ['Checkpoints already exist.  Enable resume or choose a ', ...
             'different run_name.']);
    end
end

config_path = fullfile(run_dir,'config.mat');
if exist(config_path,'file')
    saved = load(config_path,'fingerprint');
    if ~isfield(saved,'fingerprint') ...
            || ~strcmp(saved.fingerprint,fingerprint)
        error('run_omega_up_all_residual_parallel:ConfigMismatch', ...
            ['The run directory contains checkpoints for a different ', ...
             'configuration.']);
    end
else
    save(config_path,'config','fingerprint','-v7');
end

records = cell(numel(cells),1);
pending = zeros(numel(cells),1);
num_pending = 0;
num_resumed = 0;
for k = 1:numel(cells)
    checkpoint = local_checkpoint_path(checkpoint_dir,cells(k).id);
    [record,ok] = local_load_checkpoint( ...
        checkpoint,fingerprint,retry_failed,cells(k),rigorous);
    if ok
        records{k} = record;
        num_resumed = num_resumed+1;
    else
        num_pending = num_pending+1;
        pending(num_pending) = k;
    end
end
pending = pending(1:num_pending);

fprintf(['Omega_up all-functional residual run: tasks=%d, resumed=%d, ', ...
    'pending=%d, Nx=%d, Ny=%d, Ny_axis=%d\n'], ...
    numel(cells),num_resumed,numel(pending),Nx,Ny,Ny_axis);
fprintf('run directory: %s\n',run_dir);

wall_clock = tic;
worker_runtime_before = local_serial_worker_report( ...
    runtime_sha256,'omega-up-residual-before-dispatch');
worker_runtime_after = local_serial_worker_report( ...
    runtime_sha256,'omega-up-residual-after-dispatch');
if ~isempty(pending)
    use_parallel = opt.workers > 0 && exist('parpool','file') == 2;
    if opt.workers > 0 && ~use_parallel
        warning('run_omega_up_all_residual_parallel:NoParallelToolbox', ...
            'Parallel tools are unavailable; continuing serially.');
    end

    pending_records = cell(numel(pending),1);
    if use_parallel
        pool = ver10_ensure_local_pool(opt.workers);
        init_future = parfevalOnAll( ...
            pool,@omega_up_all_prepare_worker,0,project_root,mode);
        wait(init_future);
        fetchOutputs(init_future);
        worker_runtime_before = ver10_assert_worker_runtime( ...
            pool,runtime_metadata,runtime_sha256, ...
            'omega-up-residual-before-dispatch');

        cells_for_workers = cells;
        params_for_workers = mesh_params;
        parfor q = 1:numel(pending)
            k = pending(q);
            pending_records{q} = local_execute_task( ...
                cells_for_workers(k),params_for_workers, ...
                checkpoint_dir,fingerprint,rigorous);
        end
        worker_runtime_after = ver10_assert_worker_runtime( ...
            pool,runtime_metadata,runtime_sha256, ...
            'omega-up-residual-after-dispatch');
    else
        for q = 1:numel(pending)
            k = pending(q);
            fprintf('task %d/%d (id=%d, %s, ix=%d, iy=%d)\n', ...
                q,numel(pending),cells(k).id,cells(k).kind, ...
                cells(k).ix,cells(k).iy);
            pending_records{q} = local_execute_task( ...
                cells(k),mesh_params,checkpoint_dir,fingerprint,rigorous);
        end
    end

    for q = 1:numel(pending)
        records{pending(q)} = pending_records{q};
    end
end
wall_seconds = toc(wall_clock);

source_provenance = ver10_assert_source_unchanged( ...
    project_root,source_provenance);
if ~isempty(spectral_atlas_hash)
    final_atlas_hash = ver10_file_sha256( ...
        char(mesh_params.spectral_atlas_path));
    if ~strcmp(final_atlas_hash,spectral_atlas_hash)
        error('run_omega_up_all_residual_parallel:AtlasChanged', ...
            'The spectral-atlas file changed during the computation.');
    end
end
runtime_after = ver10_runtime_metadata();
runtime_after_payload = jsonencode(runtime_after);
runtime_after_sha256 = local_sha256_bytes( ...
    uint8(runtime_after_payload));
if ~strcmp(runtime_payload,runtime_after_payload) ...
        || ~strcmp(runtime_sha256,runtime_after_sha256)
    error('run_omega_up_all_residual_parallel:RuntimeChanged', ...
        ['The proof runtime changed during the computation.  Checkpoints ', ...
         'from different INTLAB, Gmsh, or engine states cannot form one ', ...
         'residual certificate.']);
end
local_assert_merged_records(records,cells,rigorous);
summary = local_summarize(records,cells,rigorous, ...
    coverage_complete,domain_coverage,wall_seconds,num_resumed,scope, ...
    source_provenance);
csv_path = fullfile(run_dir,'omega_up_all_cells.csv');
mat_path = fullfile(run_dir,'omega_up_all_results.mat');
json_path = fullfile(run_dir,'omega_up_all_manifest.json');
local_write_csv(csv_path,records);
try
    save(mat_path,'records','cells','grid','config','summary', ...
        'fingerprint','-v7.3');
catch
    save(mat_path,'records','cells','grid','config','summary', ...
        'fingerprint','-v7');
end

manifest = struct();
manifest.schema = 'lowerboundsineq.omega_up_all_residual.manifest.v2';
manifest.created_utc = local_utc_timestamp();
manifest.git_commit = source_provenance.git_commit;
manifest.source = source_provenance;
manifest.runtime = runtime_metadata;
manifest.runtime_sha256 = runtime_sha256;
manifest.runtime_rechecked_after_run = true;
manifest.worker_runtime_validation = struct( ...
    'policy',config.worker_runtime_validation_policy, ...
    'before_dispatch',worker_runtime_before, ...
    'after_dispatch',worker_runtime_after, ...
    'all_computation_workers_checked',true, ...
    'runtime_sha256',runtime_sha256);
manifest.mode = mode;
manifest.rigor = local_ternary(rigorous, ...
    'certified_interval','exploratory_double');
if strcmpi(char(mesh_params.estimator_backend),'bernstein-strong')
    manifest.estimator = ...
        'signed-strong-residual-with-global-Bernstein-bubbles';
    manifest.spectral_endpoint_method = ...
        'Lehmann-Goerisch+conforming-Ritz+CR-Liu';
    manifest.Hdiv_flux_reconstruction = false;
else
    manifest.estimator = ...
        'signed-residual-resolvent-with-shifted-RT-majorant';
    manifest.spectral_endpoint_method = 'legacy-residual-flux-backend';
    manifest.Hdiv_flux_reconstruction = true;
end
manifest.explicit_eigenfunctions = false;
manifest.high_mode_truncation = false;
manifest.single_estimator_call_per_direction_cell = true;
manifest.functionals = {'J1_lower','J2_lower','JUP_upper'};
manifest.functional_scope = scope;
manifest.domain_coverage = domain_coverage;
manifest.config = config;
manifest.summary = summary;
manifest.complete_certificate = summary.complete_certificate;
manifest.resume = struct( ...
    'enabled',resume, ...
    'retry_failed',retry_failed, ...
    'records_loaded',num_resumed, ...
    'checkpoint_directory', ...
        ver10_portable_path(project_root,checkpoint_dir), ...
    'checkpoint_granularity','one MAT file per task');
manifest.files = struct( ...
    'csv',ver10_portable_path(project_root,csv_path), ...
    'mat',ver10_portable_path(project_root,mat_path), ...
    'json',ver10_portable_path(project_root,json_path), ...
    'config_mat',ver10_portable_path(project_root,config_path));
manifest.hashes = struct( ...
    'algorithm','SHA-256', ...
    'csv',ver10_file_sha256(csv_path), ...
    'mat',ver10_file_sha256(mat_path), ...
    'config_mat',ver10_file_sha256(config_path));
if rigorous
    manifest.warning = '';
else
    manifest.warning = ...
        'Double output is exploratory and is not a proof certificate.';
end
local_write_json(json_path,manifest);

fprintf(['completed: applicable=%d, outside-disk=%d, errors=%d, ', ...
    'all-certified=%d\n'], ...
    summary.num_applicable,summary.num_skipped_outside_disk, ...
    summary.num_errors,summary.complete_certificate);
fprintf('min J1 lower = %.12g, min J2 lower = %.12g, ', ...
    summary.minimum_J1_lower,summary.minimum_J2_lower);
fprintf('max JUP upper = %.12g\n',summary.maximum_JUP_upper);
fprintf(['functional certificates: J1=%d, J2=%d, JUP=%d ', ...
    '(JUP scope: %s)\n'], ...
    summary.J1_complete_certificate, ...
    summary.J2_complete_certificate, ...
    summary.JUP_complete_certificate,scope.JUP_region);
fprintf('CSV: %s\nMAT: %s\nJSON: %s\n', ...
    csv_path,mat_path,json_path);
end


function record = local_execute_task( ...
    cell_def,mesh_params,checkpoint_dir,fingerprint,rigorous)
started = tic;
try
    record = verify_omega_up_all_functionals_cell(cell_def,mesh_params);
catch ME
    record = local_failure_record(cell_def,rigorous,ME,toc(started));
end
checkpoint = local_checkpoint_path(checkpoint_dir,cell_def.id);
local_save_checkpoint(checkpoint,record,fingerprint);
end


function record = local_failure_record(cell_def,rigorous,ME,elapsed)
record = struct();
record.task_id = double(cell_def.id);
record.kind = char(cell_def.kind);
record.ix = double(cell_def.ix);
record.iy = double(cell_def.iy);
record.x_lo = double(cell_def.x_lo);
record.x_hi = double(cell_def.x_hi);
record.y_lo = double(cell_def.y_lo);
record.y_hi = double(cell_def.y_hi);
record.x_mid = (record.x_lo+record.x_hi)/2;
record.y_mid = (record.y_lo+record.y_hi)/2;
record.direction = local_ternary( ...
    strcmp(record.kind,'rectangle'),'x','y');
record.applicable = true;
record.domain_applicable = true;
record.computed = false;
record.rigorous = logical(rigorous);
record.estimator = ...
    'signed-residual-resolvent-with-shifted-RT-majorant';
record.estimator_calls = 1;
record.J1_in_scope = logical(cell_def.J1_in_scope);
record.J2_in_scope = logical(cell_def.J2_in_scope);
record.JUP_in_scope = logical(cell_def.JUP_in_scope);
record.num_scoped_functionals = ...
    double(record.J1_in_scope)+double(record.J2_in_scope) ...
    +double(record.JUP_in_scope);
record.J1_sign_ok = false;
record.J2_sign_ok = false;
record.JUP_sign_ok = false;
record.all_signs_ok = false;
record.all_scoped_signs_ok = false;
record.J1_certified = false;
record.J2_certified = false;
record.JUP_certified = false;
record.all_certified = false;
record.elapsed_seconds = elapsed;
record.status = 'error';
record.error_identifier = ME.identifier;
record.error_message = ME.message;
end


function [cells,grid] = local_build_tasks(epsI,Nx,Ny,Ny_axis,scope)
x_left = I_intval('0.5');
y_top = sqrt(I_intval('3'))/2;
x_nodes = I_zeros(1,Nx+1);
y_nodes = I_zeros(1,Ny+1);
y_axis_nodes = I_zeros(1,Ny_axis+1);
for k = 0:Nx
    x_nodes(k+1) = x_left+2*epsI*k/Nx;
end
for k = 0:Ny
    y_nodes(k+1) = y_top-epsI*k/Ny;
end
for k = 0:Ny_axis
    y_axis_nodes(k+1) = y_top-epsI*k/Ny_axis;
end

template = struct( ...
    'id',0,'kind','','ix',0,'iy',0, ...
    'x_lo',0,'x_hi',0,'y_lo',0,'y_hi',0, ...
    'J1_in_scope',false,'J2_in_scope',false,'JUP_in_scope',false);
cells = repmat(template,Nx*Ny+Ny_axis,1);
id = 0;
for ix = 1:Nx
    for iy = 1:Ny
        id = id+1;
        cells(id).id = id;
        cells(id).kind = 'rectangle';
        cells(id).ix = ix;
        cells(id).iy = iy;
        cells(id).x_lo = I_inf(x_nodes(ix));
        cells(id).x_hi = I_sup(x_nodes(ix+1));
        cells(id).y_lo = I_inf(y_nodes(iy+1));
        cells(id).y_hi = I_sup(y_nodes(iy));
        [cells(id).J1_in_scope,cells(id).J2_in_scope, ...
            cells(id).JUP_in_scope] = local_cell_scope( ...
                scope,cells(id).y_lo,cells(id).y_hi);
    end
end
for iy = 1:Ny_axis
    id = id+1;
    cells(id).id = id;
    cells(id).kind = 'axis';
    cells(id).ix = 0;
    cells(id).iy = iy;
    cells(id).x_lo = I_inf(x_left);
    cells(id).x_hi = I_sup(x_left);
    cells(id).y_lo = I_inf(y_axis_nodes(iy+1));
    cells(id).y_hi = I_sup(y_axis_nodes(iy));
    [cells(id).J1_in_scope,cells(id).J2_in_scope, ...
        cells(id).JUP_in_scope] = local_cell_scope( ...
            scope,cells(id).y_lo,cells(id).y_hi);
end

grid = struct();
grid.x_nodes = x_nodes;
grid.y_nodes = y_nodes;
grid.y_axis_nodes = y_axis_nodes;
grid.Nx = Nx;
grid.Ny = Ny;
grid.Ny_axis = Ny_axis;
grid.functional_scope = scope;
end


function coverage = local_domain_coverage(cells,grid,scope)
% Record the region actually covered by estimator cells.  These are the
% outward cell endpoints, not merely the requested eps/threshold labels.
canonical_eps = I_intval('0.122');
canonical_x_left = I_intval('0.5');
canonical_x_right = canonical_x_left+2*canonical_eps;
canonical_y_top = sqrt(I_intval('3'))/2;
canonical_y_bottom = canonical_y_top-canonical_eps;

coverage = struct();
coverage.schema = ...
    'lowerboundsineq.omega-up-residual-domain-coverage.v1';
coverage.grid_x_left_lower = I_inf(grid.x_nodes(1));
coverage.grid_x_left_upper = I_sup(grid.x_nodes(1));
coverage.grid_x_right_lower = I_inf(grid.x_nodes(end));
coverage.grid_x_right_upper = I_sup(grid.x_nodes(end));
coverage.grid_y_bottom_lower = I_inf(grid.y_nodes(end));
coverage.grid_y_bottom_upper = I_sup(grid.y_nodes(end));
coverage.grid_y_top_lower = I_inf(grid.y_nodes(1));
coverage.grid_y_top_upper = I_sup(grid.y_nodes(1));
coverage.axis_y_bottom_lower = I_inf(grid.y_axis_nodes(end));
coverage.axis_y_bottom_upper = I_sup(grid.y_axis_nodes(end));
coverage.axis_y_top_lower = I_inf(grid.y_axis_nodes(1));
coverage.axis_y_top_upper = I_sup(grid.y_axis_nodes(1));
coverage.canonical_eps_up = I_mid(canonical_eps);
coverage.canonical_x_right_lower = I_inf(canonical_x_right);
coverage.canonical_x_right_upper = I_sup(canonical_x_right);
coverage.canonical_y_bottom_lower = I_inf(canonical_y_bottom);
coverage.canonical_y_bottom_upper = I_sup(canonical_y_bottom);

coverage.canonical_ver10_omega_up_covered = ...
    coverage.grid_x_left_lower <= I_inf(canonical_x_left) ...
    && coverage.grid_x_right_upper >= I_sup(canonical_x_right) ...
    && coverage.grid_y_bottom_lower <= I_inf(canonical_y_bottom) ...
    && coverage.grid_y_top_upper >= I_sup(canonical_y_top) ...
    && coverage.axis_y_bottom_lower <= I_inf(canonical_y_bottom) ...
    && coverage.axis_y_top_upper >= I_sup(canonical_y_top);

jup_rect = strcmp({cells.kind},'rectangle') & [cells.JUP_in_scope];
jup_axis = strcmp({cells.kind},'axis') & [cells.JUP_in_scope];
if scope.JUP_requested && any(jup_rect) && any(jup_axis)
    rect_start = min([cells(jup_rect).y_lo]);
    axis_start = min([cells(jup_axis).y_lo]);
    % Both directional estimates are needed; hence the certified common
    % band starts at the larger of their actual lower endpoints.
    coverage.JUP_rectangle_y_start = rect_start;
    coverage.JUP_axis_y_start = axis_start;
    coverage.JUP_certified_y_start_upper = max(rect_start,axis_start);
    coverage.JUP_requested_band_covered = ...
        coverage.JUP_certified_y_start_upper ...
        <= scope.JUP_y_min_lower;
else
    coverage.JUP_rectangle_y_start = NaN;
    coverage.JUP_axis_y_start = NaN;
    coverage.JUP_certified_y_start_upper = NaN;
    coverage.JUP_requested_band_covered = ~scope.JUP_requested;
end
end


function scope = local_scope_config(name,jup_y_minI)
allowed = {'split','all','lower-only','jup-band-only'};
if ~any(strcmp(name,allowed))
    error('run_omega_up_all_residual_parallel:BadFunctionalScope', ...
        ['functional_scope must be ''split'', ''all'', ', ...
         '''lower-only'', or ''jup-band-only''.']);
end
scope = struct();
scope.name = name;
scope.JUP_y_min = I_mid(jup_y_minI);
scope.JUP_y_min_lower = I_inf(jup_y_minI);
scope.JUP_y_min_upper = I_sup(jup_y_minI);
scope.J1_requested = any(strcmp(name,{'split','all','lower-only'}));
scope.J2_requested = scope.J1_requested;
scope.JUP_requested = any(strcmp(name,{'split','all','jup-band-only'}));
scope.J1_region = local_ternary( ...
    scope.J1_requested,'full_Omega_up','out_of_scope');
scope.J2_region = local_ternary( ...
    scope.J2_requested,'full_Omega_up','out_of_scope');
if strcmp(name,'all')
    scope.JUP_region = 'full_Omega_up';
elseif scope.JUP_requested
    scope.JUP_region = sprintf('Omega_up_intersect_y_ge_%.12g', ...
        scope.JUP_y_min);
else
    scope.JUP_region = 'out_of_scope';
end
end


function [J1,J2,JUP] = local_cell_scope(scope,~,y_hi)
J1 = logical(scope.J1_requested);
J2 = logical(scope.J2_requested);
if ~scope.JUP_requested
    JUP = false;
elseif strcmp(scope.name,'all')
    JUP = true;
else
    % Every cell intersecting y>=threshold is evaluated on the whole cell.
    % The resulting harmless superset keeps the curved/split boundary
    % rigorous without introducing a second estimator call.
    JUP = y_hi >= scope.JUP_y_min_lower;
end
end


function params = local_default_mesh_parameters(overrides)
params = struct( ...
    'estimator_backend','bernstein-strong', ...
    'bernstein_degree',11, ...
    'N_LG',16, ...
    'N_rho',64, ...
    'bound_order',2, ...
    'trial_order',4, ...
    'N_trial',8, ...
    'RT_order',4, ...
    'spectral_atlas_path','', ...
    'eigenfunction_error_method','residual_flux', ...
    'lambda2_analytic_lower',true, ...
    'rt_solve_strategy','midpoint-defect');
names = fieldnames(overrides);
for k = 1:numel(names)
    params.(names{k}) = overrides.(names{k});
end
end


function atlas = local_validate_spectral_atlas( ...
    filename,expected_sha256,expected_commit,rigorous)
before = ver10_file_sha256(filename);
if ~strcmp(before,expected_sha256)
    error('run_omega_up_all_residual_parallel:AtlasHashMismatch', ...
        'The spectral-atlas digest changed before loading.');
end
payload = load(filename,'atlas');
after = ver10_file_sha256(filename);
if ~strcmp(after,expected_sha256)
    error('run_omega_up_all_residual_parallel:AtlasHashMismatch', ...
        'The spectral-atlas file changed while it was loaded.');
end
if ~isfield(payload,'atlas')
    error('run_omega_up_all_residual_parallel:BadSpectralAtlas', ...
        'The spectral-atlas MAT file has no ATLAS variable.');
end
atlas = payload.atlas;
required = {'schema','rigorous','source','runtime','runtime_sha256', ...
    'runtime_rechecked_after_run','parallel_workers', ...
    'worker_runtime_validation'};
for k = 1:numel(required)
    if ~isfield(atlas,required{k})
        error('run_omega_up_all_residual_parallel:BadSpectralAtlas', ...
            'Missing atlas.%s.',required{k});
    end
end
if ~strcmp(char(atlas.schema), ...
        'lowerboundsineq.omega-up-spectral-cr-lg-atlas.v1')
    error('run_omega_up_all_residual_parallel:BadSpectralAtlas', ...
        'Unsupported spectral-atlas schema.');
end
flag = atlas.rigorous;
valid_flag = islogical(flag) && isscalar(flag);
if ~valid_flag || (rigorous && ~flag)
    error('run_omega_up_all_residual_parallel:BadSpectralAtlas', ...
        'The spectral-atlas rigor flag is invalid.');
end
if ~isstruct(atlas.source) || ~isscalar(atlas.source) ...
        || ~isfield(atlas.source,'git_commit') ...
        || ~isfield(atlas.source,'source_unchanged_after_run') ...
        || ~islogical(atlas.source.source_unchanged_after_run) ...
        || ~isscalar(atlas.source.source_unchanged_after_run) ...
        || ~strcmp(char(atlas.source.git_commit),char(expected_commit)) ...
        || ~atlas.source.source_unchanged_after_run
    error('run_omega_up_all_residual_parallel:AtlasSourceMismatch', ...
        ['The spectral atlas was not generated from the same unchanged ', ...
         'source commit as this run.']);
end
if ~isstruct(atlas.runtime) || ~isscalar(atlas.runtime) ...
        || ~(ischar(atlas.runtime_sha256) ...
            || (isstring(atlas.runtime_sha256) ...
                && isscalar(atlas.runtime_sha256))) ...
        || ~islogical(atlas.runtime_rechecked_after_run) ...
        || ~isscalar(atlas.runtime_rechecked_after_run) ...
        || ~atlas.runtime_rechecked_after_run
    error('run_omega_up_all_residual_parallel:BadAtlasRuntime', ...
        'The spectral atlas omits a rechecked full runtime snapshot.');
end
atlas_runtime_sha256 = local_sha256_bytes( ...
    uint8(jsonencode(atlas.runtime)));
if ~strcmp(char(atlas.runtime_sha256),atlas_runtime_sha256)
    error('run_omega_up_all_residual_parallel:BadAtlasRuntime', ...
        'The spectral-atlas runtime SHA-256 does not match its runtime.');
end
atlas_workers = atlas.parallel_workers;
if ~(isnumeric(atlas_workers) && isreal(atlas_workers) ...
        && isscalar(atlas_workers) && isfinite(atlas_workers) ...
        && atlas_workers >= 0 && atlas_workers == floor(atlas_workers))
    error('run_omega_up_all_residual_parallel:BadAtlasWorkers', ...
        'The spectral-atlas worker count is invalid.');
end
ver10_validate_worker_runtime_record( ...
    atlas.worker_runtime_validation,atlas_runtime_sha256, ...
    double(atlas_workers),'spectral atlas worker runtime validation');
end


function fingerprint = local_config_fingerprint(config)
m = config.mesh_params;
r = config.runtime;
fingerprint = sprintf([ ...
    'impl=%s|mode=%s|eps=[%.17g,%.17g]|dx=[%.17g,%.17g]|', ...
    'dy=[%.17g,%.17g]|Nx=%d|Ny=%d|Nya=%d|selected=%d|', ...
    'selection=%s|scope=%s|jupy=[%.17g,%.17g]|', ...
    'backend=%s|bp=%d|NLG=%d|Nrho=%d|bo=%d|to=%d|Nt=%d|RT=%d|', ...
    'ef=%s|l2=%s|rt=%s|atlas=%s|', ...
    'engine=%s|engine-version=%s|release=%s|computer=%s|', ...
    'startintlab=%s|intval=%s|intlab-tree=%s|gmsh-version=%s|gmsh=%s|', ...
    'runtime=%s|worker-runtime-policy=%s|', ...
    'git=%s|dirty=%d'], ...
    config.estimator_implementation, ...
    config.mode,config.eps_up_lower,config.eps_up_upper, ...
    config.dx_requested_lower,config.dx_requested_upper, ...
    config.dy_requested_lower,config.dy_requested_upper, ...
    config.Nx,config.Ny,config.Ny_axis, ...
    config.selected_task_count,config.selection_signature, ...
    config.functional_scope.name, ...
    config.functional_scope.JUP_y_min_lower, ...
    config.functional_scope.JUP_y_min_upper, ...
    char(m.estimator_backend),m.bernstein_degree, ...
    m.N_LG,m.N_rho,m.bound_order, ...
    m.trial_order,m.N_trial,m.RT_order, ...
    char(m.eigenfunction_error_method), ...
    local_value_string(m.lambda2_analytic_lower), ...
    char(m.rt_solve_strategy), ...
    config.spectral_atlas_sha256, ...
    char(r.engine),char(r.engine_version),char(r.release), ...
    char(r.computer),char(r.startintlab_sha256), ...
    char(r.intval_constructor_sha256),char(r.intlab_tree_sha256), ...
    char(r.gmsh_version), ...
    char(r.gmsh_binary_sha256), ...
    config.runtime_sha256, ...
    config.worker_runtime_validation_policy, ...
    config.source_git_commit,config.source_git_dirty);
end


function name = local_default_run_name(config)
m = config.mesh_params;
raw = sprintf([ ...
    'omega_up_all_%s_%s_eps%.6g_dx%.6g_dy%.6g_Nya%d_', ...
    '%s_p%d'], ...
    config.mode,config.functional_scope.name, ...
    config.eps_up,config.dx_actual,config.dy_actual,config.Ny_axis, ...
    char(m.estimator_backend),m.bernstein_degree);
name = local_safe_name(raw);
end


function name = local_safe_name(value)
name = regexprep(char(value),'[^A-Za-z0-9_.-]+','_');
if isempty(name)
    error('run_omega_up_all_residual_parallel:BadRunName', ...
        'run_name must contain at least one safe character.');
end
end


function checkpoint = local_checkpoint_path(directory,task_id)
checkpoint = fullfile(directory,sprintf('task_%06d.mat',task_id));
end


function [record,ok] = local_load_checkpoint( ...
    checkpoint,fingerprint,retry_failed,cell_def,expected_rigorous)
record = [];
ok = false;
if ~exist(checkpoint,'file')
    return;
end
try
    payload = load(checkpoint,'record','fingerprint');
    if ~isfield(payload,'record') || ~isfield(payload,'fingerprint') ...
            || ~strcmp(payload.fingerprint,fingerprint)
        return;
    end
    [coherent,successful] = local_validate_record( ...
        payload.record,cell_def,expected_rigorous);
    if ~coherent || (retry_failed && ~successful)
        return;
    end
    record = payload.record;
    ok = true;
catch
    % A checkpoint interrupted during save is recomputed.
    record = [];
    ok = false;
end
end


function local_assert_merged_records(records,cells,expected_rigorous)
if ~iscell(records) || numel(records) ~= numel(cells)
    error('run_omega_up_all_residual_parallel:BadMergedRecordCount', ...
        'The merged record count does not match the selected task count.');
end
for k = 1:numel(cells)
    [coherent,~,reason] = local_validate_record( ...
        records{k},cells(k),expected_rigorous);
    if ~coherent
        error('run_omega_up_all_residual_parallel:InvalidMergedRecord', ...
            'Task %d has an invalid merged record: %s', ...
            cells(k).id,reason);
    end
end
end


function [coherent,successful,reason] = local_validate_record( ...
    record,cell_def,expected_rigorous)
coherent = false;
successful = false;
reason = 'record is not a scalar structure';
if ~isstruct(record) || ~isscalar(record)
    return;
end

expected_kind = char(cell_def.kind);
expected_scope = logical([ ...
    cell_def.J1_in_scope,cell_def.J2_in_scope,cell_def.JUP_in_scope]);
expected_scoped_count = sum(double(expected_scope));
identity_ok = ...
    local_number_field_equals(record,'task_id',cell_def.id) ...
    && local_string_field_equals(record,'kind',expected_kind) ...
    && local_number_field_equals(record,'ix',cell_def.ix) ...
    && local_number_field_equals(record,'iy',cell_def.iy) ...
    && local_number_field_equals(record,'x_lo',cell_def.x_lo) ...
    && local_number_field_equals(record,'x_hi',cell_def.x_hi) ...
    && local_number_field_equals(record,'y_lo',cell_def.y_lo) ...
    && local_number_field_equals(record,'y_hi',cell_def.y_hi) ...
    && local_logical_field_equals( ...
        record,'J1_in_scope',expected_scope(1)) ...
    && local_logical_field_equals( ...
        record,'J2_in_scope',expected_scope(2)) ...
    && local_logical_field_equals( ...
        record,'JUP_in_scope',expected_scope(3)) ...
    && local_number_field_equals( ...
        record,'num_scoped_functionals',expected_scoped_count) ...
    && local_logical_field_equals( ...
        record,'rigorous',logical(expected_rigorous));
if ~identity_ok
    reason = 'task identity, geometry, scope, or rigor does not match';
    return;
end

if expected_scoped_count == 0
    [coherent,reason] = local_validate_skip_record( ...
        record,'skipped_out_of_scope',true,expected_rigorous);
    successful = coherent;
    return;
end
if strcmp(expected_kind,'rectangle') ...
        && local_cell_is_certified_outside_disk(cell_def)
    [coherent,reason] = local_validate_skip_record( ...
        record,'skipped_outside_disk',false,expected_rigorous);
    successful = coherent;
    return;
end

status = local_get(record,'status','');
if local_string_value_equals(status,'error')
    coherent = ...
        local_logical_field_equals(record,'applicable',true) ...
        && local_logical_field_equals(record,'domain_applicable',true) ...
        && local_logical_field_equals(record,'computed',false) ...
        && local_number_field_equals(record,'estimator_calls',1) ...
        && local_all_false(record,{ ...
            'J1_sign_ok','J2_sign_ok','JUP_sign_ok', ...
            'all_signs_ok','all_scoped_signs_ok', ...
            'J1_certified','J2_certified','JUP_certified', ...
            'all_certified'}) ...
        && local_nonnegative_finite_field(record,'elapsed_seconds');
    successful = false;
    if coherent
        reason = '';
    else
        reason = 'error record is internally inconsistent';
    end
    return;
end

[coherent,successful,reason] = local_validate_computed_record( ...
    record,expected_kind,expected_scope,expected_rigorous);
end


function [coherent,reason] = local_validate_skip_record( ...
    record,expected_status,expected_domain_applicable,expected_rigorous)
nan_fields = { ...
    'lambda_lower','lambda_upper','lambda2_lower', ...
    'dlambda_lower','dlambda_upper', ...
    'ddlambda_lower','ddlambda_upper', ...
    'J1_lower','J2_lower','JUP_upper', ...
    'R_upper','eps_D_upper','eps_a_upper','eps_0_upper', ...
    'residual_gap_lower','eigenfunction_strong_residual_upper', ...
    'material_strong_residual_upper','qmin', ...
    'point_flux_upper','gradient_transport_upper', ...
    'source_transport_upper'};
coherent = ...
    local_string_field_equals(record,'status',expected_status) ...
    && local_logical_field_equals(record,'applicable',false) ...
    && local_logical_field_equals( ...
        record,'domain_applicable',expected_domain_applicable) ...
    && local_logical_field_equals(record,'computed',false) ...
    && local_number_field_equals(record,'estimator_calls',0) ...
    && local_logical_field_equals(record,'all_signs_ok',true) ...
    && local_logical_field_equals(record,'all_scoped_signs_ok',true) ...
    && local_all_false(record,{ ...
        'J1_sign_ok','J2_sign_ok','JUP_sign_ok', ...
        'J1_certified','J2_certified','JUP_certified'}) ...
    && local_logical_field_equals( ...
        record,'all_certified',logical(expected_rigorous)) ...
    && local_fields_are_nan(record,nan_fields) ...
    && local_nonnegative_finite_field(record,'elapsed_seconds');
if coherent
    reason = '';
else
    reason = 'geometric/scope skip record is internally inconsistent';
end
end


function [coherent,successful,reason] = local_validate_computed_record( ...
    record,expected_kind,expected_scope,expected_rigorous)
coherent = false;
successful = false;
reason = 'computed record metadata is inconsistent';
expected_direction = local_ternary( ...
    strcmp(expected_kind,'rectangle'),'x','y');
if ~(local_logical_field_equals(record,'applicable',true) ...
        && local_logical_field_equals(record,'domain_applicable',true) ...
        && local_logical_field_equals(record,'computed',true) ...
        && local_number_field_equals(record,'estimator_calls',1) ...
        && local_string_field_one_of(record,'estimator',{ ...
            'signed-strong-residual-with-global-Bernstein-bubbles', ...
            'signed-residual-resolvent-with-shifted-RT-majorant'}) ...
        && local_string_field_equals(record,'direction',expected_direction) ...
        && local_nonnegative_finite_field(record,'elapsed_seconds'))
    return;
end

estimator = char(record.estimator);
if strcmp(estimator, ...
        'signed-strong-residual-with-global-Bernstein-bubbles')
    strong_fields = { ...
        'residual_gap_lower', ...
        'eigenfunction_strong_residual_upper', ...
        'material_strong_residual_upper'};
    if ~local_fields_are_finite_scalars(record,strong_fields) ...
            || ~(double(record.residual_gap_lower) > 0) ...
            || ~(double(record.eigenfunction_strong_residual_upper) >= 0) ...
            || ~(double(record.material_strong_residual_upper) >= 0)
        reason = 'Bernstein strong-residual diagnostics are invalid';
        return;
    end
end

finite_fields = { ...
    'lambda_lower','lambda_upper','lambda2_lower', ...
    'dlambda_lower','dlambda_upper', ...
    'ddlambda_upper','JUP_upper', ...
    'R_upper','eps_D_upper','eps_a_upper','eps_0_upper','qmin'};
if ~local_fields_are_finite_scalars(record,finite_fields)
    reason = 'computed spectral, residual, or functional bound is nonfinite';
    return;
end
lower_requested = expected_scope(1) || expected_scope(2);
if lower_requested
    if ~local_fields_are_finite_scalars( ...
            record,{'ddlambda_lower','J1_lower','J2_lower'})
        reason = 'an in-scope lower functional bound is nonfinite';
        return;
    end
elseif ~(local_negative_infinity_field(record,'ddlambda_lower') ...
        && local_negative_infinity_field(record,'J1_lower') ...
        && local_negative_infinity_field(record,'J2_lower'))
    reason = 'upper-only record does not carry the expected -Inf sentinels';
    return;
end

lambda_lower = double(record.lambda_lower);
lambda_upper = double(record.lambda_upper);
lambda2_lower = double(record.lambda2_lower);
dlambda_lower = double(record.dlambda_lower);
dlambda_upper = double(record.dlambda_upper);
ddlambda_lower = double(record.ddlambda_lower);
ddlambda_upper = double(record.ddlambda_upper);
if ~(lambda_lower > 0 && lambda_lower <= lambda_upper ...
        && lambda2_lower > lambda_upper ...
        && dlambda_lower <= dlambda_upper ...
        && ddlambda_lower <= ddlambda_upper ...
        && double(record.R_upper) >= 0 ...
        && double(record.eps_D_upper) >= 0 ...
        && double(record.eps_a_upper) >= 0 ...
        && double(record.eps_0_upper) >= 0 ...
        && double(record.qmin) > 0)
    reason = 'computed bound ordering or residual nonnegativity failed';
    return;
end

signs = logical([ ...
    double(record.J1_lower) > 0, ...
    double(record.J2_lower) > 0, ...
    double(record.JUP_upper) < 0]);
expected_all_signs = all(~expected_scope | signs);
expected_certified = logical(expected_rigorous) ...
    & expected_scope & signs;
expected_all_certified = logical(expected_rigorous) ...
    && expected_all_signs;
if ~(local_logical_field_equals(record,'J1_sign_ok',signs(1)) ...
        && local_logical_field_equals(record,'J2_sign_ok',signs(2)) ...
        && local_logical_field_equals(record,'JUP_sign_ok',signs(3)) ...
        && local_logical_field_equals( ...
            record,'all_signs_ok',expected_all_signs) ...
        && local_logical_field_equals( ...
            record,'all_scoped_signs_ok',expected_all_signs) ...
        && local_logical_field_equals( ...
            record,'J1_certified',logical(expected_certified(1))) ...
        && local_logical_field_equals( ...
            record,'J2_certified',logical(expected_certified(2))) ...
        && local_logical_field_equals( ...
            record,'JUP_certified',logical(expected_certified(3))) ...
        && local_logical_field_equals( ...
            record,'all_certified',expected_all_certified))
    reason = 'sign and certification flags do not match finite bounds';
    return;
end

if expected_all_certified
    expected_status = 'verified';
elseif ~expected_rigorous && expected_all_signs
    expected_status = 'exploratory_signs_hold';
else
    expected_status = 'unverified_sign';
end
if ~local_string_field_equals(record,'status',expected_status)
    reason = 'status does not match the recomputed sign/certificate state';
    return;
end

coherent = true;
successful = expected_all_signs ...
    && (~expected_rigorous || expected_all_certified);
reason = '';
end


function tf = local_cell_is_certified_outside_disk(cell_def)
radius_sq_lower = I_inf( ...
    I_intval(cell_def.x_lo)^2+I_intval(cell_def.y_lo)^2);
tf = radius_sq_lower > 1;
end


function tf = local_number_field_equals(record,name,expected)
value = local_get(record,name,[]);
tf = (isnumeric(value) || islogical(value)) ...
    && isscalar(value) && isfinite(double(value)) ...
    && double(value) == double(expected);
end


function tf = local_logical_field_equals(record,name,expected)
value = local_get(record,name,[]);
tf = (isnumeric(value) || islogical(value)) ...
    && isscalar(value) && isfinite(double(value)) ...
    && any(double(value) == [0,1]) ...
    && logical(value) == logical(expected);
end


function tf = local_string_field_equals(record,name,expected)
tf = local_string_value_equals(local_get(record,name,''),expected);
end


function tf = local_string_field_one_of(record,name,expected)
tf = false;
for k = 1:numel(expected)
    if local_string_field_equals(record,name,expected{k})
        tf = true;
        return;
    end
end
end


function tf = local_string_value_equals(value,expected)
tf = (ischar(value) || (isstring(value) && isscalar(value))) ...
    && strcmp(char(value),char(expected));
end


function tf = local_fields_are_finite_scalars(record,names)
tf = true;
for k = 1:numel(names)
    value = local_get(record,names{k},[]);
    if ~((isnumeric(value) || islogical(value)) ...
            && isscalar(value) && isfinite(double(value)))
        tf = false;
        return;
    end
end
end


function tf = local_fields_are_nan(record,names)
tf = true;
for k = 1:numel(names)
    value = local_get(record,names{k},[]);
    if ~(isnumeric(value) && isscalar(value) && isnan(double(value)))
        tf = false;
        return;
    end
end
end


function tf = local_negative_infinity_field(record,name)
value = local_get(record,name,[]);
tf = isnumeric(value) && isscalar(value) ...
    && isinf(double(value)) && double(value) < 0;
end


function tf = local_all_false(record,names)
tf = true;
for k = 1:numel(names)
    if ~local_logical_field_equals(record,names{k},false)
        tf = false;
        return;
    end
end
end


function tf = local_nonnegative_finite_field(record,name)
value = local_get(record,name,[]);
tf = isnumeric(value) && isscalar(value) ...
    && isfinite(double(value)) && double(value) >= 0;
end


function local_save_checkpoint(checkpoint,record,fingerprint)
directory = fileparts(checkpoint);
temporary = [tempname(directory),'.mat'];
cleanup = onCleanup(@() local_delete_if_present(temporary));
save(temporary,'record','fingerprint','-v7');
[ok,message] = movefile(temporary,checkpoint,'f');
if ~ok
    error('run_omega_up_all_residual_parallel:CheckpointWriteFailed', ...
        'Could not publish checkpoint %s: %s',checkpoint,message);
end
clear cleanup
end


function local_delete_if_present(filename)
if exist(filename,'file')
    delete(filename);
end
end


function summary = local_summarize( ...
    records,cells,rigorous,coverage_complete,domain_coverage,wall_seconds, ...
    num_resumed,scope,source_provenance)
n = numel(records);
record_ids = NaN(n,1);
applicable = false(n,1);
domain_applicable = false(n,1);
computed = false(n,1);
skipped = false(n,1);
skipped_scope = false(n,1);
errors = false(n,1);
all_certified = false(n,1);
J1_scope = false(n,1);
J2_scope = false(n,1);
JUP_scope = false(n,1);
J1_certified = false(n,1);
J2_certified = false(n,1);
JUP_certified = false(n,1);
for k = 1:n
    r = records{k};
    record_ids(k) = local_get(r,'task_id',NaN);
    applicable(k) = logical(local_get(r,'applicable',true));
    domain_applicable(k) = logical( ...
        local_get(r,'domain_applicable',applicable(k)));
    computed(k) = logical(local_get(r,'computed',false));
    status = local_get(r,'status','');
    skipped(k) = strcmp(status,'skipped_outside_disk');
    skipped_scope(k) = strcmp(status,'skipped_out_of_scope');
    errors(k) = strcmp(status,'error');
    all_certified(k) = logical(local_get(r,'all_certified',false));
    J1_scope(k) = logical(local_get(r,'J1_in_scope',true));
    J2_scope(k) = logical(local_get(r,'J2_in_scope',true));
    JUP_scope(k) = logical(local_get(r,'JUP_in_scope',true));
    J1_certified(k) = logical(local_get(r,'J1_certified',false));
    J2_certified(k) = logical(local_get(r,'J2_certified',false));
    JUP_certified(k) = logical(local_get(r,'JUP_certified',false));
end

expected_ids = double([cells.id].');
missing_ids = setdiff(expected_ids,record_ids);
unexpected_ids = setdiff(record_ids,expected_ids);
num_duplicate_ids = n-numel(unique(record_ids));
exact_selected_ids = isempty(missing_ids) && isempty(unexpected_ids) ...
    && num_duplicate_ids == 0 && n == numel(expected_ids);
exact_full_ids = logical(coverage_complete) && exact_selected_ids ...
    && isequal(sort(record_ids(:)),(1:n).');
source_clean = source_provenance.git_status_available ...
    && ~source_provenance.git_dirty ...
    && ~isempty(source_provenance.git_commit) ...
    && logical(local_get( ...
        source_provenance,'source_unchanged_after_run',false));
certificate_context_ok = exact_selected_ids && source_clean;
canonical_context_ok = certificate_context_ok ...
    && logical(domain_coverage.canonical_ver10_omega_up_covered);

rect = strcmp({cells.kind}.','rectangle');
axis = strcmp({cells.kind}.','axis');
neutral_or_certified = all_certified | skipped | skipped_scope;
J1_required = domain_applicable & J1_scope;
J2_required = domain_applicable & J2_scope;
JUP_required = domain_applicable & JUP_scope;
J1_requested = logical(scope.J1_requested);
J2_requested = logical(scope.J2_requested);
JUP_requested = logical(scope.JUP_requested);

J1_summary = local_function_summary( ...
    J1_requested,J1_required,J1_certified,errors, ...
    rigorous,coverage_complete,canonical_context_ok);
J2_summary = local_function_summary( ...
    J2_requested,J2_required,J2_certified,errors, ...
    rigorous,coverage_complete,canonical_context_ok);
JUP_summary = local_function_summary( ...
    JUP_requested,JUP_required,JUP_certified,errors, ...
    rigorous,coverage_complete,canonical_context_ok ...
        && logical(domain_coverage.JUP_requested_band_covered));

summary = struct();
summary.num_tasks = n;
summary.num_rectangles = sum(rect);
summary.num_axis_intervals = sum(axis);
summary.num_applicable = sum(applicable);
summary.num_computed = sum(computed);
summary.num_skipped_outside_disk = sum(skipped);
summary.num_skipped_out_of_scope = sum(skipped_scope);
summary.num_errors = sum(errors);
summary.num_resumed = num_resumed;
summary.num_unique_task_ids = numel(unique(record_ids));
summary.num_missing_task_ids = numel(missing_ids);
summary.num_unexpected_task_ids = numel(unexpected_ids);
summary.num_duplicate_task_ids = num_duplicate_ids;
summary.exact_selected_task_id_coverage = exact_selected_ids;
summary.exact_full_task_id_coverage = exact_full_ids;
summary.source_clean = source_clean;
summary.source_unchanged_after_run = logical(local_get( ...
    source_provenance,'source_unchanged_after_run',false));
summary.num_J1_required = sum(J1_required);
summary.num_J2_required = sum(J2_required);
summary.num_JUP_required = sum(JUP_required);
summary.num_J1_certified = J1_summary.num_certified;
summary.num_J2_certified = J2_summary.num_certified;
summary.num_JUP_certified = JUP_summary.num_certified;
summary.rectangles_all_certified = ...
    all(neutral_or_certified(rect));
summary.axis_all_certified = all(neutral_or_certified(axis));
summary.coverage_complete = logical(coverage_complete);
summary.canonical_ver10_omega_up_covered = logical( ...
    domain_coverage.canonical_ver10_omega_up_covered);
summary.JUP_requested_band_covered = logical( ...
    domain_coverage.JUP_requested_band_covered);
summary.JUP_certified_y_start_upper = ...
    domain_coverage.JUP_certified_y_start_upper;
summary.rigorous = logical(rigorous);
summary.J1_complete_certificate = J1_summary.complete_certificate;
summary.J2_complete_certificate = J2_summary.complete_certificate;
summary.JUP_complete_certificate = JUP_summary.complete_certificate;
summary.functionals = struct( ...
    'J1',J1_summary,'J2',J2_summary,'JUP',JUP_summary);
requested = [J1_requested,J2_requested,JUP_requested];
completed = [J1_summary.complete_certificate, ...
    J2_summary.complete_certificate,JUP_summary.complete_certificate];
summary.complete_certificate = any(requested) ...
    && canonical_context_ok ...
    && exact_full_ids ...
    && all(completed(requested));
summary.minimum_J1_lower = ...
    local_record_extreme(records,'J1_lower','min',J1_required);
summary.minimum_J2_lower = ...
    local_record_extreme(records,'J2_lower','min',J2_required);
summary.maximum_JUP_upper = ...
    local_record_extreme(records,'JUP_upper','max',JUP_required);
summary.minimum_J1_rectangle_lower = local_record_extreme( ...
    records,'J1_lower','min',J1_required & rect);
summary.minimum_J1_axis_lower = local_record_extreme( ...
    records,'J1_lower','min',J1_required & axis);
summary.minimum_J2_rectangle_lower = local_record_extreme( ...
    records,'J2_lower','min',J2_required & rect);
summary.minimum_J2_axis_lower = local_record_extreme( ...
    records,'J2_lower','min',J2_required & axis);
summary.maximum_JUP_rectangle_upper = local_record_extreme( ...
    records,'JUP_upper','max',JUP_required & rect);
summary.maximum_JUP_axis_upper = local_record_extreme( ...
    records,'JUP_upper','max',JUP_required & axis);
summary.wall_seconds = wall_seconds;
summary.sum_task_seconds = local_record_sum(records,'elapsed_seconds');
end


function result = local_function_summary( ...
    requested,required,certified,errors,rigorous,coverage_complete, ...
    certificate_context_ok)
result = struct();
result.requested = logical(requested);
result.num_required = sum(required);
result.num_certified = sum(certified & required);
result.num_failed = sum(required & (~certified | errors));
result.complete_certificate = logical(requested) ...
    && logical(rigorous) ...
    && logical(coverage_complete) ...
    && logical(certificate_context_ok) ...
    && result.num_required > 0 ...
    && result.num_certified == result.num_required ...
    && ~any(errors & required);
end


function value = local_record_extreme(records,field,kind,mask)
values = [];
for k = 1:numel(records)
    if ~mask(k)
        continue;
    end
    candidate = local_get(records{k},field,NaN);
    if isnumeric(candidate) && isscalar(candidate) && isfinite(candidate)
        values(end+1,1) = candidate; %#ok<AGROW>
    end
end
if isempty(values)
    value = NaN;
elseif strcmp(kind,'min')
    value = min(values);
else
    value = max(values);
end
end


function value = local_record_sum(records,field)
value = 0;
for k = 1:numel(records)
    candidate = local_get(records{k},field,0);
    if isnumeric(candidate) && isscalar(candidate) && isfinite(candidate)
        value = value+candidate;
    end
end
end


function local_write_csv(filename,records)
fid = fopen(filename,'w');
if fid < 0
    error('run_omega_up_all_residual_parallel:CsvOpenFailed', ...
        'Cannot open %s.',filename);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,[ ...
    'task_id,kind,ix,iy,x_lo,x_hi,y_lo,y_hi,applicable,', ...
    'domain_applicable,computed,rigorous,estimator_calls,', ...
    'J1_in_scope,J2_in_scope,JUP_in_scope,num_scoped_functionals,', ...
    'lambda_lower,lambda_upper,lambda2_lower,', ...
    'dlambda_lower,dlambda_upper,ddlambda_lower,ddlambda_upper,', ...
    'J1_lower,J1_sign_ok,J1_certified,J2_lower,J2_sign_ok,J2_certified,', ...
        'JUP_upper,JUP_sign_ok,JUP_certified,all_signs_ok,', ...
        'all_scoped_signs_ok,all_certified,', ...
        'R_upper,eps_D_upper,eps_a_upper,eps_0_upper,', ...
        'residual_gap_lower,eigenfunction_strong_residual_upper,', ...
        'material_strong_residual_upper,qmin,', ...
        'point_flux_upper,gradient_transport_upper,source_transport_upper,', ...
    'elapsed_seconds,status,error_identifier,error_message\n']);
for k = 1:numel(records)
    r = records{k};
    fprintf(fid,[ ...
        '%d,%s,%d,%d,%.17g,%.17g,%.17g,%.17g,', ...
        '%d,%d,%d,%d,%d,%d,%d,%d,%d,', ...
        '%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,', ...
        '%.17g,%d,%d,%.17g,%d,%d,%.17g,%d,%d,%d,%d,%d,', ...
        '%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,', ...
        '%.17g,%.17g,%.17g,%.17g,', ...
        '%.17g,%s,%s,%s\n'], ...
        local_get(r,'task_id',NaN), ...
        local_csv_escape(local_get(r,'kind','')), ...
        local_get(r,'ix',0),local_get(r,'iy',0), ...
        local_get(r,'x_lo',NaN),local_get(r,'x_hi',NaN), ...
        local_get(r,'y_lo',NaN),local_get(r,'y_hi',NaN), ...
        local_get(r,'applicable',false), ...
        local_get(r,'domain_applicable',false), ...
        local_get(r,'computed',false), ...
        local_get(r,'rigorous',false), ...
        local_get(r,'estimator_calls',0), ...
        local_get(r,'J1_in_scope',false), ...
        local_get(r,'J2_in_scope',false), ...
        local_get(r,'JUP_in_scope',false), ...
        local_get(r,'num_scoped_functionals',0), ...
        local_get(r,'lambda_lower',NaN), ...
        local_get(r,'lambda_upper',NaN), ...
        local_get(r,'lambda2_lower',NaN), ...
        local_get(r,'dlambda_lower',NaN), ...
        local_get(r,'dlambda_upper',NaN), ...
        local_get(r,'ddlambda_lower',NaN), ...
        local_get(r,'ddlambda_upper',NaN), ...
        local_get(r,'J1_lower',NaN), ...
        local_get(r,'J1_sign_ok',false), ...
        local_get(r,'J1_certified',false), ...
        local_get(r,'J2_lower',NaN), ...
        local_get(r,'J2_sign_ok',false), ...
        local_get(r,'J2_certified',false), ...
        local_get(r,'JUP_upper',NaN), ...
        local_get(r,'JUP_sign_ok',false), ...
        local_get(r,'JUP_certified',false), ...
        local_get(r,'all_signs_ok',false), ...
        local_get(r,'all_scoped_signs_ok',false), ...
        local_get(r,'all_certified',false), ...
        local_get(r,'R_upper',NaN), ...
        local_get(r,'eps_D_upper',NaN), ...
        local_get(r,'eps_a_upper',NaN), ...
        local_get(r,'eps_0_upper',NaN), ...
        local_get(r,'residual_gap_lower',NaN), ...
        local_get(r,'eigenfunction_strong_residual_upper',NaN), ...
        local_get(r,'material_strong_residual_upper',NaN), ...
        local_get(r,'qmin',NaN), ...
        local_get(r,'point_flux_upper',NaN), ...
        local_get(r,'gradient_transport_upper',NaN), ...
        local_get(r,'source_transport_upper',NaN), ...
        local_get(r,'elapsed_seconds',NaN), ...
        local_csv_escape(local_get(r,'status','')), ...
        local_csv_escape(local_get(r,'error_identifier','')), ...
        local_csv_escape(local_get(r,'error_message','')));
end
end


function value = local_get(record,name,default)
if isstruct(record) && isfield(record,name)
    value = record.(name);
else
    value = default;
end
end


function escaped = local_csv_escape(value)
if ischar(value)
    value = value;
elseif isstring(value)
    value = char(value);
elseif isnumeric(value) || islogical(value)
    value = num2str(value);
else
    value = char(value);
end
value = strrep(value,sprintf('\r'),' ');
value = strrep(value,sprintf('\n'),' ');
value = strrep(value,'"','""');
escaped = ['"',value,'"'];
end


function local_write_json(filename,manifest)
try
    encoded = jsonencode(manifest,'PrettyPrint',true);
catch
    encoded = jsonencode(manifest);
end
fid = fopen(filename,'w');
if fid < 0
    error('run_omega_up_all_residual_parallel:JsonOpenFailed', ...
        'Cannot open %s.',filename);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',encoded);
end


function digest = local_sha256_bytes(bytes)
try
    engine = javaMethod( ...
        'getInstance','java.security.MessageDigest','SHA-256');
catch ME
    error('run_omega_up_all_residual_parallel:NoSHA256', ...
        'Cannot initialize SHA-256: %s',ME.message);
end
bytes = uint8(bytes(:));
if ~isempty(bytes)
    engine.update(typecast(bytes,'int8'));
end
raw = typecast(engine.digest(),'uint8');
digest = lower(reshape(dec2hex(raw,2).',1,[]));
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


function stamp = local_utc_timestamp()
if exist('OCTAVE_VERSION','builtin')
    stamp = strftime('%Y-%m-%dT%H:%M:%SZ',gmtime(time()));
else
    now_utc = datetime('now','TimeZone','UTC');
    stamp = char(string(now_utc,'yyyy-MM-dd''T''HH:mm:ss''Z'''));
end
end


function value = local_decimal_value(input)
if ischar(input) || isstring(input)
    value = I_intval(char(input));
else
    value = I_intval(input);
end
end


function value = local_value_string(input)
if islogical(input)
    value = local_ternary(input,'true','false');
elseif ischar(input) || isstring(input)
    value = char(input);
else
    value = sprintf('%.17g',double(input));
end
end


function result = local_ternary(condition,yes_value,no_value)
if condition
    result = yes_value;
else
    result = no_value;
end
end


function tf = local_positive_scalar_input(x)
if ischar(x) || isstring(x)
    x = str2double(x);
end
tf = isnumeric(x) && isscalar(x) && isfinite(x) && x > 0;
end


function tf = local_scope_name_input(x)
if ~(ischar(x) || (isstring(x) && isscalar(x)))
    tf = false;
    return;
end
tf = any(strcmp(lower(char(x)), ...
    {'split','all','lower-only','jup-band-only'}));
end


function tf = local_optional_positive_integer(x)
tf = isempty(x) || (isnumeric(x) && isscalar(x) ...
    && isfinite(x) && x >= 1 && x == floor(x));
end


function tf = local_optional_task_ids(x)
tf = isempty(x) || (isnumeric(x) && isvector(x) ...
    && all(isfinite(x(:))) && all(x(:) >= 1) ...
    && all(x(:) == floor(x(:))));
end


function tf = local_nonnegative_integer(x)
tf = isnumeric(x) && isscalar(x) && isfinite(x) ...
    && x >= 0 && x == floor(x);
end
