function manifest = run_upper_conjecture_atlas(varargin)
%RUN_UPPER_CONJECTURE_ATLAS Parallel verifier for Siudeja's upper bound.
%
% Examples:
%   % Fast diagnostic (never labelled certified):
%   run_upper_conjecture_atlas('mode','double','backend','p1', ...
%       'max_cells',8,'p1_n',12);
%
%   % Rigorous run (INTLAB_ROOT must be set on every worker):
%   run_upper_conjecture_atlas('mode','interval','workers',20, ...
%       'backend','p4','trial_order',4,'N_trial',8, ...
%       'dx',5e-4,'dy',5e-4);
%
% The production atlas reaches y=0.851 and overlaps the independent
% Omega_up residual certificate, whose default lower interface is y=0.85.
% Both endpoints are configurable and recorded in the manifest.

project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(project_root, 'src', 'upper_conjecture'));
source_provenance = ver10_source_provenance(project_root);

p = inputParser;
addParameter(p, 'mode', 'double', @(s) ischar(s) || isstring(s));
addParameter(p, 'backend', 'p4', @(s) ischar(s) || isstring(s));
addParameter(p, 'p1_n', 24, @(x) isnumeric(x) && isscalar(x) && x >= 2);
addParameter(p, 'trial_order', 4, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'N_trial', 8, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 2);
addParameter(p, 'y_min', 0.06, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'y_up', 0.851, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
addParameter(p, 'residual_y_start', 0.85, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
addParameter(p, 'dx', 0.001, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'dy', 0.001, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'workers', 0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'max_cells', Inf, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'output_dir', fullfile(project_root, 'results', 'upper_conjecture'), ...
    @(s) ischar(s) || isstring(s));
parse(p, varargin{:});
opt = p.Results;
mode = lower(char(opt.mode));
backend = lower(char(opt.backend));
if ~ismember(mode, {'double','interval'})
    error('mode must be "double" or "interval".');
end
if ~ismember(backend, {'p1','p4','high_order'})
    error('backend must be "p1", "p4", or "high_order".');
end
if opt.y_min >= opt.y_up
    error('y_min must be below y_up.');
end
if opt.y_min > 0.06
    error('y_min cannot exceed the certified thin-sector endpoint 0.06.');
end
if opt.residual_y_start >= opt.y_up
    error('residual_y_start must be below y_up to create a proof overlap.');
end

upper_prepare_runtime(project_root, mode);
[cells,cover_complete] = upper_build_compact_cells( ...
    opt.y_min,opt.y_up,opt.dx,opt.dy,floor(opt.max_cells));
if strcmp(mode,'interval') && cover_complete ...
        && (~source_provenance.git_status_available ...
            || source_provenance.git_dirty ...
            || isempty(source_provenance.git_commit))
    error('run_upper_conjecture_atlas:DirtyProductionTree', ...
        ['A complete interval atlas must start from a clean Git tree. ', ...
         'Commit the implementation or use max_cells for a smoke run.']);
end
if strcmp(backend,'p1')
    mesh = upper_make_reference_mesh(opt.p1_n);
    matrices = upper_reference_matrices(mesh);
else
    mesh = [];
    matrices = upper_make_high_order_model( ...
        opt.trial_order,opt.N_trial,mode);
end

if ~isfolder(opt.output_dir)
    mkdir(opt.output_dir);
end
stamp = datestr(now, 'yyyymmdd_HHMMSS');
csv_path = fullfile(opt.output_dir, ['upper_atlas_', mode, '_', stamp, '.csv']);
json_path = fullfile(opt.output_dir, ['upper_manifest_', mode, '_', stamp, '.json']);

if opt.workers > 0 && exist('parpool', 'file') == 2
    pool = gcp('nocreate');
    if isempty(pool)
        parpool(opt.workers);
    elseif pool.NumWorkers ~= opt.workers
        delete(pool);
        parpool(opt.workers);
    end
    pool = gcp('nocreate');
    init_future = parfevalOnAll( ...
        pool,@upper_prepare_runtime,0,project_root,mode);
    wait(init_future);
    for f = 1:numel(init_future)
        if ~isempty(init_future(f).Error)
            rethrow(init_future(f).Error);
        end
    end
end

results = cell(numel(cells), 1);
use_parallel = opt.workers > 0 && exist('parpool', 'file') == 2;
if use_parallel
    parfor k = 1:numel(cells)
        upper_prepare_runtime(project_root, mode);
        results{k} = upper_verify_trial_cell(cells(k), mesh, matrices, mode);
    end
else
    for k = 1:numel(cells)
        results{k} = upper_verify_trial_cell(cells(k), mesh, matrices, mode);
    end
end

thin = upper_verify_thin_sector(mode, opt.y_min);
source_provenance = ver10_assert_source_unchanged( ...
    project_root,source_provenance);
upper_write_results_csv(csv_path, results);
csv_sha256 = ver10_file_sha256(csv_path);
margins = cellfun(@(r) r.margin_lower, results);
verified_flags = cellfun(@(r) r.verified, results);

manifest.schema = 'lowerboundsineq.upper_conjecture.v1';
if exist('OCTAVE_VERSION', 'builtin')
    manifest.created_utc = strftime('%Y-%m-%dT%H:%M:%SZ', gmtime(time()));
else
    utc_now = datetime('now', 'TimeZone', 'UTC');
    manifest.created_utc = char(string(utc_now, 'yyyy-MM-dd''T''HH:mm:ss''Z'''));
end
manifest.mode = mode;
manifest.rigor = ternary(strcmp(mode,'interval'), ...
    'certified_interval', 'exploratory_double');
manifest.source = source_provenance;
manifest.runtime = ver10_runtime_metadata();
manifest.theorem = ['lambda1 <= pi^2 L^2/(12 A^2) ', ...
    '+ sqrt(3) pi^2/(3 A)'];
manifest.normalization = 'T(x,y)=conv{(0,0),(1,0),(x,y)}, A=y/2';
manifest.thin = thin;
manifest.atlas.y_min = opt.y_min;
manifest.atlas.y_up = opt.y_up;
manifest.atlas.x_min = 0.5;
manifest.atlas.curved_boundary = 'rectangular_superset_cover';
manifest.atlas.dx = opt.dx;
manifest.atlas.dy = opt.dy;
manifest.atlas.backend = backend;
manifest.atlas.mesh_order = ternary(strcmp(backend,'p1'),1,opt.trial_order);
manifest.atlas.trial_resolution = ternary( ...
    strcmp(backend,'p1'),opt.p1_n,opt.N_trial);
manifest.atlas.number_of_cells = numel(cells);
manifest.atlas.cover_generation_complete = cover_complete;
manifest.atlas.minimum_margin_lower = min(margins);
manifest.atlas.number_verified = sum(verified_flags);
manifest.atlas.all_verified = strcmp(mode,'interval') ...
    && cover_complete && all(verified_flags);
manifest.component_complete_certificate = ...
    thin.verified && manifest.atlas.all_verified;
manifest.component_certificate_scope = sprintf( ...
    'normalized triangles with 0 < y <= %.17g',opt.y_up);
manifest.connection.residual_y_start = opt.residual_y_start;
manifest.connection.atlas_y_end = opt.y_up;
manifest.connection.overlap_y = [opt.residual_y_start,opt.y_up];
manifest.connection.required_external_certificate = ...
    'Omega_up residual Jup certificate through the equilateral point';
manifest.complete_certificate = false;
manifest.complete_certificate_reason = ...
    ['This component manifest deliberately does not assert the global ', ...
     'theorem; combine it with the overlapping residual manifest.'];
manifest.files.cells_csv = ver10_portable_path(project_root,csv_path);
manifest.files.manifest_json = ver10_portable_path(project_root,json_path);
manifest.hashes.algorithm = 'SHA-256';
manifest.hashes.cells_csv = csv_sha256;
manifest.warning = ternary(strcmp(mode,'interval'), '', ...
    'Double output is exploratory and is not a proof certificate.');

fid = fopen(json_path, 'w');
if fid < 0
    error('Cannot open manifest %s.', json_path);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', jsonencode(manifest, 'PrettyPrint', true));

fprintf('upper conjecture: mode=%s cells=%d min_margin=%.6e\n', ...
    mode, numel(cells), manifest.atlas.minimum_margin_lower);
fprintf('thin scaled margin lower: %.6e\n', thin.scaled_margin_lower);
fprintf('CSV: %s\nJSON: %s\n', csv_path, json_path);
if strcmp(mode,'interval')
    fprintf('certified cells: %d/%d; component-complete=%d\n', ...
        manifest.atlas.number_verified, numel(cells), ...
        manifest.component_complete_certificate);
else
    fprintf('EXPLORATORY DOUBLE ONLY -- not a proof certificate.\n');
end
end

function value = ternary(condition, yes_value, no_value)
if condition
    value = yes_value;
else
    value = no_value;
end
end
