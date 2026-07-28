function record = verify_omega_up_all_functionals_cell( ...
    cell_def,mesh_params)
%VERIFY_OMEGA_UP_ALL_FUNCTIONALS_CELL  One estimator call, three functionals.
%
% CELL_DEF fields:
%   id, kind ('rectangle' or 'axis'), ix, iy,
%   x_lo, x_hi, y_lo, y_hi,
% and the optional logical fields
%   J1_in_scope, J2_in_scope, JUP_in_scope.
% If the scope fields are absent, all three functionals are requested.
%
% A rectangle is tested in the x direction and an axis interval in the y
% direction.  For every applicable cell, calc_ddlami_residual_bounds is
% called exactly once with a two-sided Hessian request.  The resulting
% lambda_1 derivative enclosure is reused to certify
%
%   J1_xx,J2_xx > 0 and Jup_xx < 0,
%
% or, on x=1/2,
%
%   J1_yy,J2_yy > 0 and Jup_yy < 0.
%
% Rectangle cells certified disjoint from x^2+y^2<=1 are skipped.  A cell
% intersecting the circle is evaluated on the whole rectangle; this
% harmless superset treatment closes the curved boundary rigorously.

global INTERVAL_MODE
if isempty(INTERVAL_MODE)
    INTERVAL_MODE = 1;
end
rigorous = logical(INTERVAL_MODE);

local_validate_cell(cell_def);
[estimator_backend,bernstein_degree] = local_estimator_backend(mesh_params);
[J1_in_scope,J2_in_scope,JUP_in_scope] = ...
    local_functional_scope(cell_def);
[N_LG,N_rho,bound_order,trial_order,N_trial,RT_order] = ...
    local_mesh_parameters(mesh_params);

x_lo = local_lower_scalar(cell_def.x_lo);
x_hi = local_upper_scalar(cell_def.x_hi);
y_lo = local_lower_scalar(cell_def.y_lo);
y_hi = local_upper_scalar(cell_def.y_hi);

record = local_empty_record();
record.task_id = double(cell_def.id);
record.kind = lower(char(cell_def.kind));
record.ix = double(cell_def.ix);
record.iy = double(cell_def.iy);
record.x_lo = x_lo;
record.x_hi = x_hi;
record.y_lo = y_lo;
record.y_hi = y_hi;
record.x_mid = (x_lo+x_hi)/2;
record.y_mid = (y_lo+y_hi)/2;
record.rigorous = rigorous;
if strcmp(estimator_backend,'bernstein-strong')
    record.estimator = ...
        'signed-strong-residual-with-global-Bernstein-bubbles';
else
    record.estimator = ...
        'signed-residual-resolvent-with-shifted-RT-majorant';
end
record.J1_in_scope = J1_in_scope;
record.J2_in_scope = J2_in_scope;
record.JUP_in_scope = JUP_in_scope;
record.num_scoped_functionals = ...
    double(J1_in_scope)+double(J2_in_scope)+double(JUP_in_scope);

if record.num_scoped_functionals == 0
    record.applicable = false;
    record.domain_applicable = true;
    record.computed = false;
    record.status = 'skipped_out_of_scope';
    record.all_signs_ok = true;
    record.all_scoped_signs_ok = true;
    record.all_certified = rigorous;
    return;
end

if strcmp(record.kind,'rectangle') ...
        && local_certified_outside_disk(x_lo,y_lo)
    record.applicable = false;
    record.domain_applicable = false;
    record.computed = false;
    record.status = 'skipped_outside_disk';
    record.all_signs_ok = true;
    record.all_scoped_signs_ok = true;
    record.all_certified = rigorous;
    return;
end
record.domain_applicable = true;

xI = I_infsup(x_lo,x_hi);
yI = I_infsup(y_lo,y_hi);
if strcmp(record.kind,'rectangle')
    direction_name = 'x';
    e_direction = [I_intval(1),I_intval(0)];
else
    direction_name = 'y';
    e_direction = [I_intval(0),I_intval(1)];
end

base_triangle = local_triangle(I_mid(xI),I_mid(yI));
target_triangle = local_triangle(xI,yI);
options = local_estimator_options( ...
    mesh_params,bound_order,trial_order,N_trial);
% The lower functionals require both Hessian endpoints.  A JUP-only run may
% use the special i=1 upper path, for which the shifted material flux is not
% needed because C_minus=0.
options.upper_only = ~(J1_in_scope || J2_in_scope);

started = tic;
if strcmp(estimator_backend,'bernstein-strong')
    spectral_options = struct( ...
        'N_LG',N_LG, ...
        'N_rho',N_rho, ...
        'bound_order',bound_order);
    if isfield(mesh_params,'spectral_atlas') ...
            && ~isempty(mesh_params.spectral_atlas)
        spectral_options.reference_certificate = ...
            omega_up_spectral_atlas_certificate( ...
                mesh_params.spectral_atlas,xI,yI);
    elseif isfield(mesh_params,'spectral_atlas_path') ...
            && ~isempty(mesh_params.spectral_atlas_path)
        if ~isfield(mesh_params,'spectral_atlas_sha256') ...
                || isempty(mesh_params.spectral_atlas_sha256)
            error( ...
                'verify_omega_up_all_functionals_cell:MissingAtlasHash', ...
                'A file-backed spectral atlas requires its SHA-256 digest.');
        end
        spectral_options.reference_certificate = ...
            omega_up_spectral_atlas_certificate( ...
                mesh_params.spectral_atlas_path,xI,yI, ...
                mesh_params.spectral_atlas_sha256);
    end
    [lambda,dlambda,ddlambda_lower,ddlambda_upper,diagnostics] = ...
        calc_ddlambda1_bernstein_strong_bounds( ...
            target_triangle,e_direction,bernstein_degree, ...
            spectral_options);
else
    [lambda,dlambda,ddlambda_lower,ddlambda_upper,diagnostics] = ...
        calc_ddlami_residual_bounds( ...
            1,base_triangle,target_triangle,e_direction, ...
            N_LG,N_rho,bound_order,RT_order,options);
end
record.elapsed_seconds = toc(started);

if ~(J1_in_scope || J2_in_scope)
    % Preserve the fail-closed upper-only record contract.  The Bernstein
    % branch may compute the lower endpoint incidentally, but it is not
    % part of an upper-only certificate.
    ddlambda_lower = I_intval(-Inf);
end

curvatures = omega_up_all_functional_curvatures( ...
    direction_name,xI,yI,dlambda,ddlambda_lower,ddlambda_upper);

record.applicable = true;
record.computed = true;
record.estimator_calls = 1;
record.direction = direction_name;
record.lambda_lower = I_inf(lambda);
record.lambda_upper = I_sup(lambda);
record.dlambda_lower = I_inf(dlambda);
record.dlambda_upper = I_sup(dlambda);
record.ddlambda_lower = I_inf(ddlambda_lower);
record.ddlambda_upper = I_sup(ddlambda_upper);
record.J1_lower = I_inf(curvatures.J1.bound);
record.J2_lower = I_inf(curvatures.J2.bound);
record.JUP_upper = I_sup(curvatures.JUP.bound);
record.J1_sign_ok = logical(curvatures.J1.sign_verified);
record.J2_sign_ok = logical(curvatures.J2.sign_verified);
record.JUP_sign_ok = logical(curvatures.JUP.sign_verified);
record.all_signs_ok = ...
    (~J1_in_scope || record.J1_sign_ok) ...
    && (~J2_in_scope || record.J2_sign_ok) ...
    && (~JUP_in_scope || record.JUP_sign_ok);
record.all_scoped_signs_ok = record.all_signs_ok;
record.J1_certified = rigorous && J1_in_scope && record.J1_sign_ok;
record.J2_certified = rigorous && J2_in_scope && record.J2_sign_ok;
record.JUP_certified = rigorous && JUP_in_scope && record.JUP_sign_ok;
record.all_certified = rigorous && record.all_signs_ok;

if isfield(diagnostics.core,'R')
    record.R_upper = I_sup(diagnostics.core.R);
else
    record.R_upper = I_sup(diagnostics.core.material_residual_l2);
end
record.eps_D_upper = I_sup(diagnostics.core.eps_D);
record.eps_a_upper = I_sup(diagnostics.eps_a);
record.eps_0_upper = I_sup(diagnostics.eps_0);
record.lambda2_lower = I_inf(diagnostics.lambda_bounds(2));
record.qmin = I_inf(diagnostics.eig_factor_lower);
if isfield(diagnostics,'residual_gap')
    record.residual_gap_lower = I_inf(diagnostics.residual_gap);
end
if isfield(diagnostics,'eigenfunction_strong_residual')
    record.eigenfunction_strong_residual_upper = ...
        I_sup(diagnostics.eigenfunction_strong_residual);
end
if isfield(diagnostics,'material_strong_residual')
    record.material_strong_residual_upper = ...
        I_sup(diagnostics.material_strong_residual);
end
if isfield(diagnostics.flux,'point_residual')
    record.point_flux_upper = ...
        I_sup(diagnostics.flux.point_residual);
end
if isfield(diagnostics.flux,'gradient_transport')
    record.gradient_transport_upper = ...
        I_sup(diagnostics.flux.gradient_transport);
end
if isfield(diagnostics.flux,'source_transport')
    record.source_transport_upper = ...
        I_sup(diagnostics.flux.source_transport);
end

if record.all_certified
    record.status = 'verified';
elseif ~rigorous && record.all_signs_ok
    record.status = 'exploratory_signs_hold';
else
    record.status = 'unverified_sign';
end
end


function [backend,degree] = local_estimator_backend(params)
if isfield(params,'estimator_backend') ...
        && ~isempty(params.estimator_backend)
    backend = lower(char(params.estimator_backend));
else
    backend = 'shifted-rt';
end
switch backend
    case {'bernstein','bernstein-strong','strong-residual'}
        backend = 'bernstein-strong';
    case {'shifted-rt','rt','residual-flux'}
        backend = 'shifted-rt';
    otherwise
        error('verify_omega_up_all_functionals_cell:BadEstimatorBackend', ...
            'Unknown mesh_params.estimator_backend: %s.',backend);
end
if isfield(params,'bernstein_degree')
    degree = params.bernstein_degree;
else
    degree = 11;
end
if ~isscalar(degree) || ~isfinite(degree) ...
        || degree ~= floor(degree) || degree < 3
    error('verify_omega_up_all_functionals_cell:BadBernsteinDegree', ...
        'mesh_params.bernstein_degree must be an integer at least three.');
end
end


function options = local_estimator_options( ...
    params,bound_order,trial_order,N_trial)
options = struct( ...
    'bound_order',bound_order, ...
    'trial_order',trial_order, ...
    'N_trial',N_trial);
forward = {'eigenfunction_error_method', ...
    'lambda2_analytic_lower','rt_solve_strategy'};
for k = 1:numel(forward)
    name = forward{k};
    if isfield(params,name)
        options.(name) = params.(name);
    end
end
end


function [N_LG,N_rho,bound_order,trial_order,N_trial,RT_order] = ...
    local_mesh_parameters(params)
if ~isstruct(params) ...
        || ~isfield(params,'N_LG') || ~isfield(params,'N_rho')
    error('verify_omega_up_all_functionals_cell:BadMeshParameters', ...
        'mesh_params.N_LG and mesh_params.N_rho are required.');
end
N_LG = params.N_LG;
N_rho = params.N_rho;
if isfield(params,'bound_order')
    bound_order = params.bound_order;
elseif isfield(params,'fem_ord_LG')
    bound_order = params.fem_ord_LG;
else
    error('verify_omega_up_all_functionals_cell:BadMeshParameters', ...
        'mesh_params.bound_order is required.');
end
if isfield(params,'trial_order')
    trial_order = params.trial_order;
else
    trial_order = bound_order;
end
if isfield(params,'N_trial')
    N_trial = params.N_trial;
else
    N_trial = N_LG;
end
if isfield(params,'RT_order')
    RT_order = params.RT_order;
else
    RT_order = trial_order;
end

values = [N_LG,N_rho,bound_order,trial_order,N_trial,RT_order];
if any(~isfinite(values)) || any(values < 1) ...
        || any(values ~= floor(values))
    error('verify_omega_up_all_functionals_cell:BadMeshParameters', ...
        'All mesh resolutions and polynomial orders must be positive integers.');
end
if RT_order ~= trial_order
    error('verify_omega_up_all_functionals_cell:OrderMismatch', ...
        'The two-sided residual estimator requires RT_order=trial_order.');
end
end


function local_validate_cell(cell_def)
required = {'id','kind','ix','iy','x_lo','x_hi','y_lo','y_hi'};
for k = 1:numel(required)
    if ~isfield(cell_def,required{k})
        error('verify_omega_up_all_functionals_cell:BadCell', ...
            'Missing cell_def.%s.',required{k});
    end
end
kind = lower(char(cell_def.kind));
if ~any(strcmp(kind,{'rectangle','axis'}))
    error('verify_omega_up_all_functionals_cell:BadCell', ...
        'cell_def.kind must be ''rectangle'' or ''axis''.');
end
x_lo = local_lower_scalar(cell_def.x_lo);
x_hi = local_upper_scalar(cell_def.x_hi);
y_lo = local_lower_scalar(cell_def.y_lo);
y_hi = local_upper_scalar(cell_def.y_hi);
if ~all(isfinite([x_lo,x_hi,y_lo,y_hi])) ...
        || x_lo > x_hi || y_lo <= 0 || y_lo > y_hi ...
        || x_lo < 0.5
    error('verify_omega_up_all_functionals_cell:BadCell', ...
        'Invalid Omega_up cell endpoints.');
end
if strcmp(kind,'axis') && ~(x_lo <= 0.5 && 0.5 <= x_hi)
    error('verify_omega_up_all_functionals_cell:BadAxisCell', ...
        'An axis task must contain x=1/2.');
end
end


function [J1,J2,JUP] = local_functional_scope(cell_def)
names = {'J1_in_scope','J2_in_scope','JUP_in_scope'};
values = true(1,3);
for k = 1:numel(names)
    if isfield(cell_def,names{k})
        candidate = cell_def.(names{k});
        if ~(isscalar(candidate) ...
                && (islogical(candidate) || isnumeric(candidate)))
            error('verify_omega_up_all_functionals_cell:BadScope', ...
                'cell_def.%s must be a logical scalar.',names{k});
        end
        values(k) = logical(candidate);
    end
end
J1 = values(1);
J2 = values(2);
JUP = values(3);
end


function tf = local_certified_outside_disk(x_lo,y_lo)
% x,y are positive on Omega_up, hence x^2+y^2 is coordinatewise
% increasing.  Its lower-left value rigorously decides disjointness.
radius_sq_lower = I_inf(I_intval(x_lo)^2+I_intval(y_lo)^2);
tf = radius_sq_lower > 1;
end


function triangle = local_triangle(x,y)
triangle = I_zeros(1,6);
triangle(3) = I_intval(1);
triangle(5) = x;
triangle(6) = y;
end


function value = local_lower_scalar(x)
if isa(x,'intval')
    value = double(I_inf(x));
else
    value = double(x);
end
end


function value = local_upper_scalar(x)
if isa(x,'intval')
    value = double(I_sup(x));
else
    value = double(x);
end
end


function record = local_empty_record()
nan_value = NaN;
record = struct( ...
    'task_id',nan_value, ...
    'kind','', ...
    'ix',nan_value, ...
    'iy',nan_value, ...
    'x_lo',nan_value, ...
    'x_hi',nan_value, ...
    'y_lo',nan_value, ...
    'y_hi',nan_value, ...
    'x_mid',nan_value, ...
    'y_mid',nan_value, ...
	    'direction','', ...
	    'applicable',true, ...
	    'domain_applicable',true, ...
	    'computed',false, ...
	    'rigorous',false, ...
	    'estimator','', ...
	    'estimator_calls',0, ...
	    'J1_in_scope',true, ...
	    'J2_in_scope',true, ...
	    'JUP_in_scope',true, ...
	    'num_scoped_functionals',3, ...
    'lambda_lower',nan_value, ...
    'lambda_upper',nan_value, ...
    'lambda2_lower',nan_value, ...
    'dlambda_lower',nan_value, ...
    'dlambda_upper',nan_value, ...
    'ddlambda_lower',nan_value, ...
    'ddlambda_upper',nan_value, ...
    'J1_lower',nan_value, ...
    'J2_lower',nan_value, ...
    'JUP_upper',nan_value, ...
    'J1_sign_ok',false, ...
    'J2_sign_ok',false, ...
	    'JUP_sign_ok',false, ...
	    'all_signs_ok',false, ...
	    'all_scoped_signs_ok',false, ...
    'J1_certified',false, ...
    'J2_certified',false, ...
    'JUP_certified',false, ...
    'all_certified',false, ...
    'R_upper',nan_value, ...
    'eps_D_upper',nan_value, ...
    'eps_a_upper',nan_value, ...
    'eps_0_upper',nan_value, ...
    'residual_gap_lower',nan_value, ...
    'eigenfunction_strong_residual_upper',nan_value, ...
    'material_strong_residual_upper',nan_value, ...
    'qmin',nan_value, ...
    'point_flux_upper',nan_value, ...
    'gradient_transport_upper',nan_value, ...
    'source_transport_upper',nan_value, ...
    'elapsed_seconds',0, ...
    'status','pending', ...
    'error_identifier','', ...
    'error_message','');
end
