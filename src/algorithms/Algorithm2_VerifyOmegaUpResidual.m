function [is_verified,results,diagnostics] = ...
    Algorithm2_VerifyOmegaUpResidual( ...
        conjecture_type,eps_up,Nx,Ny,Ny_axis,mesh_params)
%ALGORITHM2_VERIFYOMEGAUPRESIDUAL  CR/LG plus strong-residual Omega_up proof.
%
% For each rectangle (x direction) and symmetry-axis interval (y
% direction), this algorithm uses exactly one degree-p Bernstein Hessian
% enclosure.  The Hessian step needs neither a high-mode cutoff nor an
% H(div) flux.  Its only spectral data are
%
%   L1_LG <= lambda_1 <= U1_Ritz < L2_CR <= lambda_2.
%
% If MESH_PARAMS.spectral_atlas_path is present, the three endpoints are
% loaded from the nearest certified anchor and transported by exact affine
% metric factors.  Otherwise CR/LG is evaluated at each cell midpoint.

global INTERVAL_MODE
if isempty(INTERVAL_MODE)
    INTERVAL_MODE = 1;
end
conjecture_type = upper(char(conjecture_type));
if ~any(strcmp(conjecture_type,{'J1','J2','JUP'}))
    error('Algorithm2_VerifyOmegaUpResidual:BadFunctional', ...
        'CONJECTURE_TYPE must be J1, J2, or JUP.');
end
if nargin < 2 || isempty(eps_up), eps_up = I_intval('0.122'); end
if nargin < 3 || isempty(Nx), Nx = 1; end
if nargin < 4 || isempty(Ny), Ny = 1; end
if nargin < 5 || isempty(Ny_axis), Ny_axis = 1; end
if nargin < 6 || ~isstruct(mesh_params)
    error('Algorithm2_VerifyOmegaUpResidual:MissingParameters', ...
        'MESH_PARAMS is required.');
end
mesh_params = local_defaults(mesh_params);
if ~isempty(mesh_params.spectral_atlas_path) ...
        && (~isfield(mesh_params,'spectral_atlas_sha256') ...
            || isempty(mesh_params.spectral_atlas_sha256))
    mesh_params.spectral_atlas_sha256 = ...
        ver10_file_sha256(char(mesh_params.spectral_atlas_path));
end
local_positive_integers([Nx,Ny,Ny_axis]);

eps_up = I_intval(eps_up);
x_nodes = I_intval('0.5')+(2*eps_up)*(0:Nx)/Nx;
y_top = sqrt(I_intval('3'))/2;
y_nodes = y_top-eps_up*(0:Ny)/Ny;
y_axis_nodes = y_top-eps_up*(0:Ny_axis)/Ny_axis;

rectangle = local_empty_report(Ny,Nx);
axis = local_empty_report(Ny_axis,1);
failures = {};
task_id = 0;
for ix = 1:Nx
    for iy = 1:Ny
        task_id = task_id+1;
        cell_def = local_cell( ...
            task_id,'rectangle',ix,iy, ...
            x_nodes(ix),x_nodes(ix+1), ...
            y_nodes(iy+1),y_nodes(iy),conjecture_type);
        [record,failure] = local_run_cell(cell_def,mesh_params);
        rectangle.records{iy,ix} = record;
        rectangle.ok(iy,ix) = local_requested_ok(record,conjecture_type);
        rectangle.bounds(iy,ix) = ...
            local_requested_bound(record,conjecture_type);
        if ~isempty(failure), failures{end+1} = failure; end %#ok<AGROW>
    end
end
for iy = 1:Ny_axis
    task_id = task_id+1;
    cell_def = local_cell( ...
        task_id,'axis',1,iy,I_intval('0.5'),I_intval('0.5'), ...
        y_axis_nodes(iy+1),y_axis_nodes(iy),conjecture_type);
    [record,failure] = local_run_cell(cell_def,mesh_params);
    axis.records{iy} = record;
    axis.ok(iy) = local_requested_ok(record,conjecture_type);
    axis.bounds(iy) = local_requested_bound(record,conjecture_type);
    if ~isempty(failure), failures{end+1} = failure; end %#ok<AGROW>
end

is_verified = logical(INTERVAL_MODE) ...
    && all(rectangle.ok(:)) && all(axis.ok(:));
rectangle.x_nodes = x_nodes;
rectangle.y_nodes = y_nodes;
axis.y_nodes = y_axis_nodes;
rectangle.extremal_bound = ...
    local_extremum(rectangle.bounds,conjecture_type);
axis.extremal_bound = local_extremum(axis.bounds,conjecture_type);

results = struct( ...
    'region','Omega_up', ...
    'functional',conjecture_type, ...
    'rectangle',rectangle, ...
    'axis',axis, ...
    'failures',{failures}, ...
    'is_verified',is_verified);
diagnostics = struct( ...
    'estimator','Bernstein-strong-residual-with-CR-LG-endpoints', ...
    'rigorous',logical(INTERVAL_MODE), ...
    'explicit_eigenfunctions',false, ...
    'high_mode_cutoff',false, ...
    'Hdiv_flux_reconstruction',false, ...
    'spectral_atlas_used', ...
        ~isempty(mesh_params.spectral_atlas_path) ...
        || (isfield(mesh_params,'spectral_atlas') ...
            && ~isempty(mesh_params.spectral_atlas)), ...
    'num_failures',numel(failures), ...
    'Nx',Nx,'Ny',Ny,'Ny_axis',Ny_axis);
end


function params = local_defaults(params)
defaults = struct( ...
    'estimator_backend','bernstein-strong', ...
    'bernstein_degree',11, ...
    'N_LG',16, ...
    'N_rho',64, ...
    'bound_order',2, ...
    'trial_order',2, ...
    'N_trial',1, ...
    'RT_order',2, ...
    'spectral_atlas_path','');
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(params,names{k}) || isempty(params.(names{k}))
        params.(names{k}) = defaults.(names{k});
    end
end
if ~strcmpi(char(params.estimator_backend),'bernstein-strong')
    error('Algorithm2_VerifyOmegaUpResidual:BadBackend', ...
        'Algorithm 2 requires estimator_backend=bernstein-strong.');
end
end


function cell_def = local_cell( ...
    id,kind,ix,iy,x_lo,x_hi,y_lo,y_hi,functional)
cell_def = struct( ...
    'id',id,'kind',kind,'ix',ix,'iy',iy, ...
    'x_lo',I_inf(x_lo),'x_hi',I_sup(x_hi), ...
    'y_lo',I_inf(y_lo),'y_hi',I_sup(y_hi), ...
    'J1_in_scope',strcmp(functional,'J1'), ...
    'J2_in_scope',strcmp(functional,'J2'), ...
    'JUP_in_scope',strcmp(functional,'JUP'));
end


function [record,failure] = local_run_cell(cell_def,params)
failure = [];
try
    record = verify_omega_up_all_functionals_cell(cell_def,params);
catch ME
    record = struct( ...
        'status','error','applicable',true, ...
        'J1_sign_ok',false,'J2_sign_ok',false,'JUP_sign_ok',false, ...
        'J1_lower',NaN,'J2_lower',NaN,'JUP_upper',NaN);
    failure = struct( ...
        'task_id',cell_def.id, ...
        'kind',cell_def.kind,'ix',cell_def.ix,'iy',cell_def.iy, ...
        'identifier',ME.identifier,'message',ME.message);
end
end


function tf = local_requested_ok(record,functional)
if isfield(record,'applicable') && ~record.applicable
    tf = true;
    return;
end
switch functional
    case 'J1', tf = logical(record.J1_sign_ok);
    case 'J2', tf = logical(record.J2_sign_ok);
    case 'JUP', tf = logical(record.JUP_sign_ok);
end
end


function value = local_requested_bound(record,functional)
if isfield(record,'applicable') && ~record.applicable
    if strcmp(functional,'JUP'), value = -Inf; else, value = Inf; end
    return;
end
switch functional
    case 'J1', value = record.J1_lower;
    case 'J2', value = record.J2_lower;
    case 'JUP', value = record.JUP_upper;
end
end


function report = local_empty_report(nrow,ncol)
report = struct( ...
    'bounds',NaN(nrow,ncol), ...
    'ok',false(nrow,ncol), ...
    'records',{cell(nrow,ncol)});
end


function value = local_extremum(bounds,functional)
if strcmp(functional,'JUP')
    value = max(bounds(:));
else
    value = min(bounds(:));
end
end


function local_positive_integers(values)
if any(~isfinite(values)) || any(values < 1) ...
        || any(values ~= floor(values))
    error('Algorithm2_VerifyOmegaUpResidual:BadGrid', ...
        'Nx, Ny, and Ny_axis must be positive integers.');
end
end
