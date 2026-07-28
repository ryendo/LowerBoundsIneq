function [lami, dlami, ddlam_lower, ddlam_upper, diagnostics] = ...
    calc_ddlami_residual_bounds(i, base_triangle, triangle, e_direction, ...
    N_LG, N_rho, fem_ord_LG, RT_order, options)
%CALC_DDLAMI_RESIDUAL_BOUNDS  Shifted-flux bounds for lambda_i''.
%
% This estimator supports every *certifiably simple* Dirichlet eigenvalue.
% It does not truncate an eigenfunction expansion.  Instead it
%
%   1. encloses lambda_{i-1}, lambda_i, lambda_{i+1} and the eigenspace
%      error with the existing Lehmann--Goerisch/Liu--Vejchodsky code;
%   2. solves one bordered conforming FE problem for a material derivative;
%   3. equilibrates its source with one shifted Raviart--Thomas solve; and
%   4. applies the energy-dual resolvent constants C_minus and C_plus.
%
% INPUT INTERFACE
%   The BASE_TRIANGLE/TRIANGLE split is retained for drop-in use by the
%   Omega_up caller.  The actual certification is assembled on TRIANGLE,
%   not merely at BASE_TRIANGLE.  If TRIANGLE has nonzero interval width, a
%   verified point eigensolve and shifted-RT residual at its midpoint are
%   transported by explicit affine metric factors.  The remaining scalar
%   forms are evaluated directly on the affine interval cell.  Hence a
%   successful INTLAB run is uniform over that target cell.
%
% OUTPUTS
%   lami         enclosure of lambda_i
%   dlami        enclosure of the first directional derivative
%   ddlam_lower  certified scalar lower bound for lambda_i''
%   ddlam_upper  certified scalar upper bound for lambda_i''
%   diagnostics  all discrete terms, error terms, gaps and flux data
%
% The formula implemented by residual_hessian_enclosure is
%
%   c_h-2J_h-eps_D-2*C_plus*R^2 <= lambda_i''
%     <= c_h-2J_h+eps_D+2*C_minus*R^2,
%
% where C_minus=0 for i=1.
%
% OPTIONS
%   bound_order  polynomial degree used by the verified eigenvalue bounds
%   trial_order  polynomial degree used for u_h and the bordered v_h solve
%   N_trial      point-mesh resolution for that trial space
%   upper_only   for i=1, return only the upper bound and skip the RT solve
%   eigenfunction_error_method
%                'residual_flux' (default for i=1) certifies the reference
%                trial eigenfunction by an equilibrated RT residual;
%                'eigenvalue_width' retains the Liu--Vejchodsky/Rayleigh
%                conversion used previously
%   lambda2_analytic_lower
%                true (default for i=1) combines the Krahn--Szego bound
%                with the Lu--Rowlett triangle fundamental-gap bound;
%                false disables both; a positive scalar supplies another
%                certified lower bound valid on target and reference
%   rt_solve_strategy
%                'midpoint-defect' (default) or 'verified-schur'
%
% Keeping BOUND_ORDER and TRIAL_ORDER independent is useful in verified
% runs: the robust P2 lower-bound pipeline can certify the adjacent gap,
% while a pointwise P3 trial sharply reduces the Hessian discretization
% error.  RT_ORDER must equal TRIAL_ORDER whenever both bounds are needed.

global INTERVAL_MODE
if nargin < 9 || isempty(options)
    options = struct();
end
upper_only = isfield(options,'upper_only') ...
    && logical(options.upper_only);
bound_order = fem_ord_LG;
if isfield(options,'bound_order') && ~isempty(options.bound_order)
    bound_order = options.bound_order;
end
trial_order = fem_ord_LG;
if isfield(options,'trial_order') && ~isempty(options.trial_order)
    trial_order = options.trial_order;
end
N_trial = N_LG;
if isfield(options,'N_trial') && ~isempty(options.N_trial)
    N_trial = options.N_trial;
end
if isfield(options,'eigenfunction_error_method') ...
        && ~isempty(options.eigenfunction_error_method)
    eigenfunction_error_method = lower(char( ...
        options.eigenfunction_error_method));
elseif i == 1
    eigenfunction_error_method = 'residual_flux';
else
    eigenfunction_error_method = 'eigenvalue_width';
end
if strcmp(eigenfunction_error_method,'subspace')
    eigenfunction_error_method = 'eigenvalue_width';
end
if ~any(strcmp(eigenfunction_error_method, ...
        {'residual_flux','eigenvalue_width'}))
    error('calc_ddlami_residual_bounds:BadEigenfunctionErrorMethod', ...
        ['options.eigenfunction_error_method must be ', ...
         '''residual_flux'' or ''eigenvalue_width''.']);
end
if strcmp(eigenfunction_error_method,'residual_flux') && i ~= 1
    error('calc_ddlami_residual_bounds:ResidualFluxNeedsFirstEigenvalue', ...
        ['The one-sided residual-flux eigenspace certificate currently ', ...
         'applies only to i=1.']);
end
if isfield(options,'lambda2_analytic_lower')
    lambda2_analytic_lower_option = ...
        options.lambda2_analytic_lower;
else
    lambda2_analytic_lower_option = (i == 1);
end
if isfield(options,'rt_solve_strategy') ...
        && ~isempty(options.rt_solve_strategy)
    rt_solve_strategy = lower(char(options.rt_solve_strategy));
else
    rt_solve_strategy = 'midpoint-defect';
end
if nargin < 8 || isempty(RT_order)
    RT_order = trial_order;
end
if upper_only && i ~= 1
    error('calc_ddlami_residual_bounds:UpperOnlyNeedsFirstEigenvalue', ...
        'The flux-free upper-only path is valid only for i=1.');
end
if ~upper_only && RT_order ~= trial_order
    error('calc_ddlami_residual_bounds:OrderMismatch', ...
        'Exact source equilibration requires RT_order=trial_order.');
end
if i < 1 || i ~= floor(i)
    error('calc_ddlami_residual_bounds:BadIndex', ...
        'i must be a positive integer.');
end

base_triangle = I_intval(base_triangle);
triangle = I_intval(triangle);
e_direction = I_intval(e_direction);
if length(base_triangle) ~= 6 || length(triangle) ~= 6
    error('calc_ddlami_residual_bounds:BadTriangle', ...
        'Triangles must be [x1 y1 x2 y2 x3 y3].');
end
canonical_prefix = [0,0,1,0];
for k = 1:4
    if I_inf(base_triangle(k)) ~= canonical_prefix(k) ...
            || I_sup(base_triangle(k)) ~= canonical_prefix(k) ...
            || I_inf(triangle(k)) ~= canonical_prefix(k) ...
            || I_sup(triangle(k)) ~= canonical_prefix(k)
        error('calc_ddlami_residual_bounds:NoncanonicalTriangle', ...
            ['This implementation requires triangles ', ...
             'conv{(0,0),(1,0),(x,y)}.']);
    end
end
if length(e_direction) ~= 2
    error('calc_ddlami_residual_bounds:BadDirection', ...
        'e_direction must have two components.');
end
is_uniform = local_has_interval_width(triangle);

% One adjacent enclosure above lambda_i is indispensable and sufficient;
% no high-mode cutoff enters the estimator.
num_eigs = i + 1;
residual_first_path = ...
    (i == 1) && strcmp(eigenfunction_error_method,'residual_flux');
if residual_first_path
    % Certifying lambda_1 by LG needs only a lower shift below lambda_2.
    % Asking LG for lambda_2 would additionally require separation from
    % lambda_3, which fails at the equilateral double eigenvalue and is
    % irrelevant to this estimator.
    num_verified_eigs = 1;
else
    num_verified_eigs = num_eigs;
end

% Do the verified eigenvalue solve only on a point triangle.  Passing a
% nonzero-width affine cell directly to the CR verified eigensolver can make
% its rough-lower-bound stage fail even for modest cell widths.  True
% eigenvalue enclosures are transported to the target cell below.
reference_triangle = I_mid(triangle);
if residual_first_path
    [lams_reference, uh_bound, K_bound, M_bound, ...
        A_xx_bound, A_xy_bound, A_yy_bound, ...
        uh_full_bound, meshCG_bound] = ...
        calc_eigen_bounds_any_order_1k_wh( ...
            num_verified_eigs,reference_triangle, ...
            N_LG,N_rho,bound_order);
    verified_eigenvalue_diagnostics = struct( ...
        'method','Lehmann-Goerisch-first-eigenvalue-only', ...
        'requires_upper_cluster_dimension',false);
else
    [lams_reference,uh_bound,K_bound,M_bound, ...
        A_xx_bound,A_xy_bound,A_yy_bound, ...
        uh_full_bound,meshCG_bound, ...
        verified_eigenvalue_diagnostics] = ...
        calc_eigen_bounds_cluster_robust( ...
            num_verified_eigs,reference_triangle, ...
            N_LG,N_rho,bound_order);
end
% Older copies of calc_eigen_bounds_any_order_1k_wh transpose lamLow before
% I_hull(lamCG,lamLow).  With implicit expansion this produces an n-by-n
% matrix whose diagonal contains the intended enclosures.
if size(lams_reference,1) == num_verified_eigs ...
        && size(lams_reference,2) == num_verified_eigs
    lams_reference = diag(lams_reference);
end
lams_reference = lams_reference(:);
if length(lams_reference) ~= num_verified_eigs
    error('calc_ddlami_residual_bounds:BadEigenvalueShape', ...
        'The existing eigenvalue enclosure routine returned an unexpected shape.');
end

% The verified eigenvalue enclosure and the conforming trial space are
% deliberately independent.  A P2 verified LG enclosure can therefore be
% combined with a sharper P3 trial without asking the lower-bound routine
% to certify P3 eigenvalues.
if ~residual_first_path ...
        && trial_order == bound_order && N_trial == N_LG
    uh_list = uh_bound;
    K_reference = K_bound;
    M_reference = M_bound;
    A_xx_reference = A_xx_bound;
    A_xy_reference = A_xy_bound;
    A_yy_reference = A_yy_bound;
    uh_full_list = uh_full_bound;
    meshCG_reference = meshCG_bound;
else
    if trial_order == bound_order && N_trial == N_LG
        meshCG_reference = meshCG_bound;
    else
        meshCG_reference = make_mesh_by_gmsh( ...
            I_mid(reference_triangle(5)),I_mid(reference_triangle(6)), ...
            1/N_trial);
    end
    [~,uh_list,uh_full_list,K_reference,M_reference, ...
        A_xx_reference,A_xy_reference,A_yy_reference,~] = ...
        laplace_eig_lagrange_detailed( ...
            trial_order,meshCG_reference.nodes,meshCG_reference.edges, ...
            meshCG_reference.elements,meshCG_reference.boundary_edges, ...
            num_eigs);
end

lams_h_reference = I_intval(zeros(num_eigs,1));
for k = 1:num_eigs
    denom_reference = uh_list(:,k)'*M_reference*uh_list(:,k);
    lams_h_reference(k) = ...
        (uh_list(:,k)'*K_reference*uh_list(:,k))/denom_reference;
end
trial_ritz_reference = verified_ritz_enclosures( ...
    uh_list,K_reference,M_reference,num_eigs);

% Uniform affine transport of the exact eigenvalues.  If S maps the
% reference triangle to the target triangle, the pullback metric is
% S^{-1}S^{-T}.
x_reference = reference_triangle(5);
y_reference = reference_triangle(6);
x_target = triangle(5);
y_target = triangle(6);
oneI = I_intval(1);
zeroI = I_intval(0);
S_inv = [oneI,(x_reference-x_target)/y_target; ...
         zeroI,y_reference/y_target];
detS = y_target/y_reference;
[eig_factor_lower,eig_factor_upper] = ...
    local_affine_metric_eigenvalue_factors(S_inv);

% Strengthen the generally conservative verified lambda_2 lower endpoint
% before affine transport.  We combine the planar Krahn--Szego inequality
% with the Lu--Rowlett fundamental gap of a triangle,
%
%   lambda_2-lambda_1 >= 64*pi^2/(9*diameter^2).
lams_reference_verified = lams_reference;
lambda2_analytic_reference = I_intval(0);
lambda2_analytic_target = I_intval(0);
lambda2_analytic_reference_info = struct('used',false);
lambda2_analytic_target_info = struct('used',false);
lambda2_upper_reference = I_intval(NaN);
if i == 1
    [lambda2_analytic_reference,lambda2_analytic_reference_info] = ...
        local_lambda2_analytic_lower( ...
            lambda2_analytic_lower_option, ...
            I_intval(x_reference),I_intval(y_reference), ...
            lams_reference(1));
    if residual_first_path
        lambda2_upper_reference = ...
            I_intval(I_sup(trial_ritz_reference(2)));
        lams_reference(2,1) = local_interval_from_lower_upper( ...
            lambda2_analytic_reference,lambda2_upper_reference, ...
            'reference lambda_2');
    else
        lams_reference(2) = local_raise_interval_lower( ...
            lams_reference(2),lambda2_analytic_reference, ...
            'reference lambda_2');
    end
end
lams = local_scale_positive_eigen_bounds( ...
    lams_reference,eig_factor_lower,eig_factor_upper);
if i == 1
    [lambda2_analytic_target,lambda2_analytic_target_info] = ...
        local_lambda2_analytic_lower( ...
            lambda2_analytic_lower_option, ...
            x_target,y_target,lams(1));
    lams(2) = local_raise_interval_lower( ...
        lams(2),lambda2_analytic_target,'target lambda_2');
end
lambda2_reference_lower = I_intval(I_inf(lams_reference(2)));
lambda2_target_lower = I_intval(I_inf(lams(2)));

% Reuse the point mesh topology and fixed coefficient vectors.  The target
% stiffness-component matrices are transported analytically from the point
% matrices through the same 2-by-2 affine map.  This is algebraically
% identical to elementwise reassembly, but it keeps every cell uncertainty
% in a handful of scalar interval coefficients instead of repeating it in
% every global matrix entry.
meshCG = local_affine_map_triangle_mesh( ...
    meshCG_reference,I_mid(x_reference),I_mid(y_reference), ...
    x_target,y_target);
[A_xx,A_xy,A_yy] = local_transport_gradient_components( ...
    A_xx_reference,A_xy_reference,A_yy_reference,S_inv,detS);
K = A_xx+A_yy;
M = detS*M_reference;

lams_h = local_scale_positive_eigen_bounds( ...
    lams_h_reference,eig_factor_lower,eig_factor_upper);

% Reference-normalized fixed trial eigenfunction.  This normalization is
% used both by the residual-flux certificate and by the affine target-cell
% normalization below.
raw_u = uh_list(:,i);
raw_u_full = uh_full_list(:,i);
mass_raw_reference = raw_u'*M_reference*raw_u;
if ~(I_inf(mass_raw_reference) > 0)
    error('calc_ddlami_residual_bounds:BadReferenceDiscreteMass', ...
        'The reference trial eigenvector has no certified positive mass.');
end
normalizer_reference = 1/sqrt(mass_raw_reference);
u_reference = raw_u*normalizer_reference;
u_full_reference = raw_u_full*normalizer_reference;
mu_reference = lams_h_reference(i);

% Only overlap, rather than a heuristic relative gap, is used here.  The
% routine refuses to label an eigenfunction unless the requested enclosure
% is a certified singleton.
clusters = auto_cluster_eigenvalues(lams, 0);
cluster_id = [];
for k = 1:numel(clusters)
    if any(clusters{k} == i)
        cluster_id = k;
        break;
    end
end
if isempty(cluster_id) || ~isequal(clusters{cluster_id}, i)
    error('calc_ddlami_residual_bounds:NotCertifiedSimple', ...
        ['lambda_%d is not a singleton in the certified enclosures. ', ...
         'Refine the eigenvalue meshes; assigning an eigenfunction would ', ...
         'not be rigorous.'], i);
end

% A direct target-cell call to calc_grad_error_bounds contains the width of
% the transported eigenvalue enclosure under a square root.  Its error is
% therefore O(sqrt(cell width)), even though the eigenspace itself changes
% by only O(cell width).  For a nonzero-width lambda_1 cell we first
% certify the point-reference eigenspace error, then transport that error
% with a two-dimensional metric perturbation estimate.  Point evaluations
% and arbitrary higher eigenvalues retain the original path.
clusters_reference = {};
reference_cluster_id = [];
delta_a_sq_reference = I_intval(NaN);
delta_b_sq_reference = I_intval(NaN);
eigfun_info_reference = struct();
eps_a_reference = I_intval(NaN);
eps_0_reference = I_intval(NaN);
transport_info = struct('used',false);
reference_eigenfunction_flux = struct('used',false);
reference_residual_certificate = struct('used',false);

if strcmp(eigenfunction_error_method,'residual_flux')
    clusters_reference = auto_cluster_eigenvalues(lams_reference,0);
    reference_cluster_id = ...
        local_find_singleton_cluster(clusters_reference,1);
    if isempty(reference_cluster_id)
        error('calc_ddlami_residual_bounds:ReferenceNotCertifiedSimple', ...
            ['lambda_1 is not a singleton at the reference point. ', ...
             'Refine the verified eigenvalue meshes.']);
    end

    % With P_e=I, v_h=0, alpha_h=mu and lambda_h=0 the existing shifted
    % RT solve gives
    %
    %   div sigma_h = -mu u_h,
    %   rho >= ||sigma_h-grad u_h||,
    %
    % which is exactly the energy-dual eigenproblem residual.
    zero_full_reference = I_zeros(size(u_full_reference,1),1);
    reference_flux_options = struct( ...
        'solve_strategy',rt_solve_strategy, ...
        'lambda1_lower',I_inf(lams_reference(1)));
    [rho_reference,reference_eigenfunction_flux] = ...
        RT_shifted_equilibrated_flux_dirichlet( ...
            meshCG_reference,trial_order,trial_order, ...
            u_full_reference,zero_full_reference,I_eye(2,2), ...
            mu_reference,I_intval(0),reference_flux_options);
    reference_eigenfunction_flux.used = true;
    [eps_a_reference,eps_0_reference, ...
        reference_residual_certificate] = ...
        local_residual_flux_first_eigenfunction_errors( ...
            lams_reference,mu_reference,rho_reference, ...
            lambda2_reference_lower);

    if is_uniform
        [eps_a,eps_0,transport_info] = ...
            local_transport_first_eigenfunction_errors( ...
                eps_a_reference,eps_0_reference, ...
                lams_reference,lams,S_inv, ...
                eig_factor_lower,eig_factor_upper, ...
                lambda2_target_lower);
        eigenfunction_error_scope = ...
            'reference-residual-flux-plus-affine-metric-transport';
    else
        eps_a = eps_a_reference;
        eps_0 = eps_0_reference;
        transport_info = struct( ...
            'used',false,'reason','point target equals reference');
        eigenfunction_error_scope = 'reference-residual-flux';
    end

    delta_a_sq_reference = ...
        I_intval(NaN(numel(clusters_reference),1));
    delta_b_sq_reference = ...
        I_intval(NaN(numel(clusters_reference),1));
    eigfun_info_reference = struct( ...
        'method','residual_flux','ok',true, ...
        'rho',rho_reference);
    delta_a_sq = delta_a_sq_reference;
    delta_b_sq = delta_b_sq_reference;
    eigfun_info = eigfun_info_reference;
    error_cluster_id = reference_cluster_id;
elseif is_uniform && i == 1
    clusters_reference = auto_cluster_eigenvalues(lams_reference,0);
    reference_cluster_id = ...
        local_find_singleton_cluster(clusters_reference,i);
    if isempty(reference_cluster_id)
        error('calc_ddlami_residual_bounds:ReferenceNotCertifiedSimple', ...
            ['lambda_1 is not a singleton at the reference point. ', ...
             'Refine the verified eigenvalue meshes.']);
    end
    [~,delta_a_sq_reference,delta_b_sq_reference, ...
        eigfun_info_reference] = calc_grad_error_bounds( ...
            lams_reference,lams_h_reference,uh_list, ...
            K_reference,M_reference,clusters_reference);
    if ~eigfun_info_reference.ok(reference_cluster_id) ...
            || isinf(I_sup(delta_b_sq_reference(reference_cluster_id)))
        error('calc_ddlami_residual_bounds:ReferenceEigenfunctionErrorFailed', ...
            'Could not certify the reference eigenspace error for lambda_1.');
    end
    [eps_a_reference,eps_0_reference] = ...
        local_phase_aligned_eigenfunction_errors( ...
            lams_reference,lams_h_reference, ...
            delta_b_sq_reference(reference_cluster_id),i);
    [eps_a,eps_0,transport_info] = ...
        local_transport_first_eigenfunction_errors( ...
            eps_a_reference,eps_0_reference,lams_reference,lams, ...
            S_inv,eig_factor_lower,eig_factor_upper, ...
            lambda2_target_lower);

    % Preserve the established diagnostics fields while explicitly marking
    % that these are reference-point subspace estimates.
    delta_a_sq = delta_a_sq_reference;
    delta_b_sq = delta_b_sq_reference;
    eigfun_info = eigfun_info_reference;
    error_cluster_id = reference_cluster_id;
    eigenfunction_error_scope = ...
        'reference-point-plus-affine-metric-transport';
else
    [~,delta_a_sq,delta_b_sq,eigfun_info] = ...
        calc_grad_error_bounds(lams,lams_h,uh_list,K,M,clusters);
    if ~eigfun_info.ok(cluster_id) ...
            || isinf(I_sup(delta_b_sq(cluster_id)))
        error('calc_ddlami_residual_bounds:EigenfunctionErrorFailed', ...
            'Could not certify the singleton eigenspace error for lambda_%d.', i);
    end
    [eps_a,eps_0] = local_phase_aligned_eigenfunction_errors( ...
        lams,lams_h,delta_b_sq(cluster_id),i);
    error_cluster_id = cluster_id;
    eigenfunction_error_scope = 'direct-target';
end

x = x_target;
y = y_target;
a = e_direction(1);
b = e_direction(2);
if ~(I_inf(y) > 0)
    error('calc_ddlami_residual_bounds:BadHeight', ...
        'The target triangle height must be certified positive.');
end

P_e = [zeroI, -a/y; -a/y, -2*b/y];
P_ee = [2*a^2/y^2, 4*a*b/y^2; ...
        4*a*b/y^2, 6*b^2/y^2];
A_e = local_shape_matrix(P_e, A_xx, A_xy, A_yy);
A_ee = local_shape_matrix(P_ee, A_xx, A_xy, A_yy);

% Work on the fixed midpoint reference triangle.  The target mass form is
% then independent of the parameter and all cell dependence is contained
% in three 2-by-2 matrices
%
%   Q=S^{-1}S^{-T}, D=S^{-1}P_eS^{-T}, H=S^{-1}P_eeS^{-T}.
%
% The physical trial functions are u_ref/sqrt(det S) and
% v_ref/sqrt(det S).  This centered representation is algebraically
% identical to physical-domain assembly, but avoids repeated occurrences
% of det(S) and hence suppresses interval dependency.
Q_metric = S_inv*S_inv';
D_metric = S_inv*P_e*S_inv';
H_metric = S_inv*P_ee*S_inv';
identity2 = I_eye(2,2);

mass_uu = u_reference'*M_reference*u_reference;
if ~(I_inf(mass_uu) > 0)
    error('calc_ddlami_residual_bounds:BadDiscreteMass', ...
        'The discrete eigenvector has no certified positive mass.');
end
sqrt_detS = sqrt(detS);
u_h = u_reference/sqrt_detS;
u_full = u_full_reference/sqrt_detS;

cov_uu = local_reference_gradient_covariance( ...
    u_reference,u_reference,A_xx_reference,A_xy_reference,A_yy_reference);
K_uu = local_covariance_metric_form(cov_uu,Q_metric);
Ae_uu = local_covariance_metric_form(cov_uu,D_metric);
Aee_uu = local_covariance_metric_form(cov_uu,H_metric);
lambda_h = K_uu/mass_uu;
alpha_h = Ae_uu/mass_uu;

% A fixed midpoint bordered solve produces v_ref.  It need not be a
% verified linear solve: v_ref is only a trial function, and the complete
% source-solve residual is bounded below by a verified RT majorant.
y_ref = I_intval(y_reference);
P_e_reference = [zeroI,-a/y_ref; ...
                 -a/y_ref,-2*b/y_ref];
A_e_reference = local_shape_matrix( ...
    P_e_reference,A_xx_reference,A_xy_reference,A_yy_reference);
K0 = I_mid(K_reference);
M0 = I_mid(M_reference);
Ae0 = I_mid(A_e_reference);
u0 = I_mid(u_reference);
u0 = u0/sqrt(u0'*M0*u0);
lambda0 = (u0'*K0*u0)/(u0'*M0*u0);
alpha0 = (u0'*Ae0*u0)/(u0'*M0*u0);
m0u = M0*u0;
bordered = [K0-lambda0*M0, m0u; m0u', 0];
rhs = [-Ae0*u0+alpha0*m0u; 0];
bordered_solution = bordered\rhs;
v0 = bordered_solution(1:end-1);
border_multiplier = bordered_solution(end);
v_reference = I_intval(v0);
v_h = v_reference/sqrt_detS;

[inside_dofs, boundary_dofs, ndof] = local_dirichlet_dofs( ...
    meshCG, trial_order);
if numel(inside_dofs) ~= numel(v0) || ndof ~= size(u_full,1)
    error('calc_ddlami_residual_bounds:DofMapMismatch', ...
        'Could not reconcile the bordered and full Lagrange DOF maps.');
end
v_full_reference = I_zeros(ndof,1);
v_full_reference(inside_dofs) = v_reference;
v_full = v_full_reference/sqrt_detS;

cov_uv = local_reference_gradient_covariance( ...
    u_reference,v_reference, ...
    A_xx_reference,A_xy_reference,A_yy_reference);
cov_vv = local_reference_gradient_covariance( ...
    v_reference,v_reference, ...
    A_xx_reference,A_xy_reference,A_yy_reference);
mass_uv = u_reference'*M_reference*v_reference;
mass_vv = v_reference'*M_reference*v_reference;
K_vv = local_covariance_metric_form(cov_vv,Q_metric);
Ae_uv = local_covariance_metric_form(cov_uv,D_metric);

scalar_forms = struct();
scalar_forms.mass_u = mass_uu;
scalar_forms.alpha_h = alpha_h;
scalar_forms.c_h = Aee_uu/mass_uu;
scalar_forms.ell_h = -Ae_uv+alpha_h*mass_uv;
scalar_forms.b_h_vv = K_vv-lambda_h*mass_vv;
scalar_forms.grad_u_sq = K_uu;
scalar_forms.grad_v_sq = K_vv;
scalar_forms.l2_v_sq = mass_vv;

if upper_only
    flux_residual = I_intval(0);
    flux_info = struct( ...
        'skipped',true, ...
        'reason','i=1 upper bound has C_minus=0');
else
    % Certify the residual only once at the midpoint, where the RT solve
    % is sharp.  For the whole affine cell, transport that dual residual
    % through the fixed-reference forms.  For
    %
    % r_p-r_0 =
    %  -a_{D-D0}(u_ref,.)-a_{Q-I}(v_ref,.)
    %  +(alpha_p-alpha_0)m(u_ref,.)
    %  +(mu_p-mu_0)m(v_ref,.),
    %
    % coercivity a_Q >= q_min a_I and Poincare give the fully explicit
    % bound below.  It is O(diam(cell)) and contains no large interval
    % mixed solve.
    mass_uu_reference = u_reference'*M_reference*u_reference;
    K_uu_reference = u_reference'*K_reference*u_reference;
    Ae_uu_reference = u_reference'*A_e_reference*u_reference;
    lambda_h_reference_point = K_uu_reference/mass_uu_reference;
    alpha_h_reference_point = Ae_uu_reference/mass_uu_reference;
    material_flux_options = struct( ...
        'solve_strategy',rt_solve_strategy, ...
        'lambda1_lower',I_inf(lams_reference(1)));
    [point_flux_residual, point_flux_info] = ...
        RT_shifted_equilibrated_flux_dirichlet( ...
            meshCG_reference,RT_order,trial_order, ...
            u_full_reference,v_full_reference,P_e_reference, ...
            alpha_h_reference_point,lambda_h_reference_point, ...
            material_flux_options);

    qmin = I_intval(I_inf(eig_factor_lower));
    if ~(I_inf(qmin) > 0)
        error('calc_ddlami_residual_bounds:BadResidualTransportMetric', ...
            'The transported residual requires q_min>0.');
    end
    Q_difference = Q_metric-identity2;
    D_difference = D_metric-P_e_reference;
    Q_difference_norm = I_intval(I_sup(norm(Q_difference,2)));
    D_difference_norm = I_intval(I_sup(norm(D_difference,2)));
    grad_u_reference = I_intval(I_sup(sqrt( ...
        I_intval(I_sup(K_uu_reference)))));
    grad_v_reference = I_intval(I_sup(sqrt( ...
        I_intval(I_sup(v_reference'*K_reference*v_reference)))));
    l2_u_reference = I_intval(I_sup(sqrt( ...
        I_intval(I_sup(mass_uu_reference)))));
    l2_v_reference = I_intval(I_sup(sqrt( ...
        I_intval(I_sup(mass_vv)))));
    alpha_difference = local_interval_abs_upper( ...
        alpha_h-alpha_h_reference_point);
    lambda_difference_discrete = local_interval_abs_upper( ...
        lambda_h-lambda_h_reference_point);
    L1_reference = I_intval(I_inf(lams_reference(1)));

    gradient_transport = ...
        D_difference_norm*grad_u_reference ...
        +Q_difference_norm*grad_v_reference;
    source_transport = ...
        (alpha_difference*l2_u_reference ...
         +lambda_difference_discrete*l2_v_reference) ...
        /sqrt(L1_reference);
    flux_residual_enclosure = ...
        (I_intval(I_sup(point_flux_residual)) ...
         +gradient_transport+source_transport)/sqrt(qmin);
    flux_residual = I_intval(I_sup(flux_residual_enclosure));
    flux_info = struct();
    flux_info.solve_strategy = 'midpoint-RT-plus-reference-transport';
    flux_info.point = point_flux_info;
    flux_info.point_residual = point_flux_residual;
    flux_info.qmin = qmin;
    flux_info.Q_difference_norm = Q_difference_norm;
    flux_info.D_difference_norm = D_difference_norm;
    flux_info.gradient_transport = gradient_transport;
    flux_info.source_transport = source_transport;
    flux_info.residual_sq_upper = flux_residual^2;
end

core = struct();
core.i = i;
core.lambda_bounds = lams;
core.lambda_h = lambda_h;
core.lambda1_lower = I_inf(lams(1));
core.u_h = u_h;
core.v_h = v_h;
core.K = K;
core.M = M;
core.A_e = A_e;
core.A_ee = A_ee;
core.P_e = P_e;
core.P_ee = P_ee;
core.eps_a = eps_a;
core.eps_0 = eps_0;
core.flux_residual = flux_residual;
core.scalar_forms = scalar_forms;
core.upper_only = upper_only;
[ddlam_lower, ddlam_upper, core_info] = ...
    residual_hessian_enclosure(core);

lami = lams(i);
dlami = I_hull(alpha_h-core_info.eps_alpha, ...
                alpha_h+core_info.eps_alpha);

if is_uniform
    enclosure_scope = 'uniform-over-target-affine-cell';
else
    enclosure_scope = 'pointwise-target-triangle';
end

diagnostics = struct();
diagnostics.estimator = ...
    'signed-residual-resolvent-with-shifted-RT-majorant';
diagnostics.rigorous = logical(INTERVAL_MODE);
diagnostics.enclosure_scope = enclosure_scope;
diagnostics.is_uniform_over_target_cell = is_uniform;
diagnostics.base_triangle = base_triangle;
diagnostics.target_triangle = triangle;
diagnostics.reference_triangle = reference_triangle;
diagnostics.reference_choice = 'target-cell-midpoint';
diagnostics.nonzero_width_interval_eigensolve_used = false;
diagnostics.point_verified_eigensolve_used = logical(INTERVAL_MODE);
diagnostics.eig_factor_lower = eig_factor_lower;
diagnostics.eig_factor_upper = eig_factor_upper;
diagnostics.lambda_bounds_reference = lams_reference;
diagnostics.lambda_bounds_reference_verified = ...
    lams_reference_verified;
diagnostics.lambda_h_reference = lams_h_reference;
diagnostics.trial_ritz_reference = trial_ritz_reference;
diagnostics.lambda2_analytic_lower_option = ...
    lambda2_analytic_lower_option;
diagnostics.lambda2_analytic_reference = ...
    lambda2_analytic_reference;
diagnostics.lambda2_analytic_target = lambda2_analytic_target;
diagnostics.lambda2_analytic_reference_info = ...
    lambda2_analytic_reference_info;
diagnostics.lambda2_analytic_target_info = ...
    lambda2_analytic_target_info;
diagnostics.lambda2_reference_lower = lambda2_reference_lower;
diagnostics.lambda2_target_lower = lambda2_target_lower;
diagnostics.eigenfunction_error_method = eigenfunction_error_method;
diagnostics.mu_reference = mu_reference;
diagnostics.reference_u_h = u_reference;
diagnostics.reference_eigenfunction_flux = ...
    reference_eigenfunction_flux;
diagnostics.reference_residual_certificate = ...
    reference_residual_certificate;
diagnostics.target_apex = [x,y];
diagnostics.e_direction = e_direction;
diagnostics.num_eigenvalues_used = num_eigs;
diagnostics.num_verified_eigenvalues_used = num_verified_eigs;
diagnostics.verified_eigenvalue_diagnostics = ...
    verified_eigenvalue_diagnostics;
diagnostics.lambda2_upper_reference = lambda2_upper_reference;
diagnostics.high_mode_truncation = false;
diagnostics.bound_order = bound_order;
diagnostics.trial_order = trial_order;
diagnostics.N_LG = N_LG;
diagnostics.N_rho = N_rho;
diagnostics.N_trial = N_trial;
diagnostics.RT_order = RT_order;
diagnostics.rt_solve_strategy = rt_solve_strategy;
diagnostics.i = i;
diagnostics.lambda_bounds = lams;
diagnostics.lambda_h_bounds = lams_h;
diagnostics.clusters = clusters;
diagnostics.singleton_cluster_id = cluster_id;
diagnostics.eigenfunction_error_scope = eigenfunction_error_scope;
diagnostics.error_cluster_id = error_cluster_id;
diagnostics.delta_a_sq = delta_a_sq(error_cluster_id);
diagnostics.delta_b_sq = delta_b_sq(error_cluster_id);
diagnostics.eigfun_info = eigfun_info;
diagnostics.clusters_reference = clusters_reference;
diagnostics.reference_cluster_id = reference_cluster_id;
diagnostics.delta_a_sq_reference = delta_a_sq_reference;
diagnostics.delta_b_sq_reference = delta_b_sq_reference;
diagnostics.eigfun_info_reference = eigfun_info_reference;
diagnostics.eps_a_reference = eps_a_reference;
diagnostics.eps_0_reference = eps_0_reference;
diagnostics.eigenfunction_transport = transport_info;
diagnostics.eps_a = eps_a;
diagnostics.eps_0 = eps_0;
diagnostics.u_h = u_h;
diagnostics.v_h = v_h;
diagnostics.boundary_dofs = boundary_dofs;
diagnostics.lambda_h = lambda_h;
diagnostics.alpha_h = alpha_h;
diagnostics.P_e = P_e;
diagnostics.P_ee = P_ee;
diagnostics.border_multiplier = border_multiplier;
diagnostics.border_residual_norm = norm( ...
    bordered*bordered_solution-rhs,2);
diagnostics.flux = flux_info;
diagnostics.core = core_info;
diagnostics.lami = lami;
diagnostics.dlami = dlami;
diagnostics.ddlam_lower = ddlam_lower;
diagnostics.ddlam_upper = ddlam_upper;
end


function cluster_id = local_find_singleton_cluster(clusters,i)
cluster_id = [];
for k = 1:numel(clusters)
    if any(clusters{k} == i)
        if isequal(clusters{k},i)
            cluster_id = k;
        end
        return;
    end
end
end


function [eps_a,eps_0] = local_phase_aligned_eigenfunction_errors( ...
    lams,lams_h,delta_b_sq,i)
% Convert a singleton subspace-angle estimate into phase-aligned H1/L2
% errors.  Only outward upper endpoints leave this routine.
db2 = min(max(I_sup(delta_b_sq),0),1);
cos_arg = I_intval(1)-I_intval(db2);
cos_arg_lower = max(I_inf(cos_arg),0);
cos_lower = I_inf(sqrt(I_intval(cos_arg_lower)));

eps_0_sq = I_intval(2)-2*I_intval(cos_lower);
eps_0 = I_intval(I_sup(sqrt( ...
    I_intval(max(I_sup(eps_0_sq),0)))));

eps_a_sq = I_intval(I_sup(lams(i))) ...
    +I_intval(I_sup(lams_h(i))) ...
    -2*I_intval(I_inf(lams(i)))*I_intval(cos_lower);
eps_a = I_intval(I_sup(sqrt( ...
    I_intval(max(I_sup(eps_a_sq),0)))));
end


function [eps_a,eps_0,info] = ...
    local_residual_flux_first_eigenfunction_errors( ...
    lams_reference,mu_reference,rho_reference,L2_reference_lower)
% Reference-point eigenfunction certificate from one equilibrated flux.
%
% For an L2-normalized conforming u_h with Rayleigh quotient mu, let
% div(sigma)=-mu*u_h and rho=||sigma-grad u_h||.  The energy-dual residual
% is at most rho.  Resolving it against the certified reference gap gives
%
%   eta = L2/(L2-U_mu) rho,
%   q^2 = eta^2/L2,
%
% followed by the phase-aligned H1/L2 errors below.

if length(lams_reference) < 2
    error('calc_ddlami_residual_bounds:ResidualFluxNeedsGap', ...
        'The residual-flux certificate requires lambda_2.');
end

L2_reference = I_intval(I_inf(L2_reference_lower));
U1_reference = I_intval(I_sup(lams_reference(1)));
mu_upper = I_intval(I_sup(mu_reference));
rho = I_intval(I_sup(rho_reference));
if ~(I_inf(L2_reference) > 0) ...
        || ~(I_inf(U1_reference) > 0) ...
        || ~(I_inf(mu_upper) > 0) ...
        || I_inf(rho) < 0 || ~isfinite(I_sup(rho))
    error('calc_ddlami_residual_bounds:BadResidualFluxData', ...
        'The residual-flux certificate received invalid spectral data.');
end

gap = L2_reference-mu_upper;
if ~(I_inf(gap) > 0)
    error('calc_ddlami_residual_bounds:ResidualFluxGapNotCertified', ...
        ['The residual-flux certificate requires ', ...
         'L_2(reference)>upper(mu).']);
end

eta_enclosure = L2_reference/gap*rho;
eta = I_intval(I_sup(eta_enclosure));
q2_enclosure = eta^2/L2_reference;
q2_upper = I_sup(q2_enclosure);
if ~isfinite(q2_upper) || q2_upper < 0 || ~(q2_upper < 1)
    error('calc_ddlami_residual_bounds:ResidualFluxAngleNotCertified', ...
        ['The residual-flux certificate requires ', ...
         'eta^2/L_2(reference)<1.']);
end
q2 = I_intval(max(q2_upper,0));

cos_phase = sqrt(I_intval(1)-q2);
cos_phase_lower = I_inf(cos_phase);
phase_loss = I_intval(max(I_sup( ...
    I_intval(1)-I_intval(cos_phase_lower)),0));

eps_a_sq = eta^2+U1_reference*phase_loss^2;
eps_a = I_intval(I_sup(sqrt( ...
    I_intval(max(I_sup(eps_a_sq),0)))));
eps_0_sq = 2*phase_loss;
eps_0 = I_intval(I_sup(sqrt( ...
    I_intval(max(I_sup(eps_0_sq),0)))));

info = struct();
info.used = true;
info.L2_reference = L2_reference;
info.U1_reference = U1_reference;
info.mu = mu_reference;
info.mu_upper = mu_upper;
info.gap = gap;
info.rho = rho;
info.eta = eta;
info.q2 = q2;
info.cos_phase_lower = cos_phase_lower;
info.eps_a = eps_a;
info.eps_0 = eps_0;
end


function [eps_a,eps_0,info] = ...
    local_transport_first_eigenfunction_errors( ...
    eps_a_reference,eps_0_reference,lams_reference,lams_target, ...
    S_inv,qmin_bound,qmax_bound,L2_target_lower)
% Transport a certified reference lambda_1 eigenspace to an affine cell.
%
% On the fixed reference triangle let
%
%   Q = S^{-1}S^{-T},  E = Q-I.
%
% The residual of the transported reference eigenfunction is controlled
% linearly in ||E||.  A Davis--Kahan/resolvent estimate then controls its
% component orthogonal to the target first eigenspace.  Unlike inserting
% cell-wide eigenvalue intervals into a Rayleigh identity, every added
% term below is O(cell width).

if length(lams_reference) < 1 || length(lams_target) < 2
    error('calc_ddlami_residual_bounds:TransportNeedsAdjacentGap', ...
        'Transport of lambda_1 requires reference lambda_1 and target lambda_2.');
end

qmin = I_intval(I_inf(qmin_bound));
qmax = I_intval(I_sup(qmax_bound));
if ~(I_inf(qmin) > 0) || ~isfinite(I_sup(qmax))
    error('calc_ddlami_residual_bounds:BadTransportMetric', ...
        'The affine metric must have finite positive eigenvalue bounds.');
end

Q = S_inv*S_inv';
E = Q-I_eye(2,2);
E_norm = I_intval(I_sup(norm(E,2)));

U1_reference = I_intval(I_sup(lams_reference(1)));
L2_target = I_intval(I_inf(L2_target_lower));
U1_target = I_intval(I_sup(lams_target(1)));
if ~(I_inf(U1_reference) > 0) || ~(I_inf(L2_target) > 0)
    error('calc_ddlami_residual_bounds:BadTransportSpectrum', ...
        'The transported eigenspace estimate requires positive spectral bounds.');
end

R0_enclosure = E_norm/sqrt(qmin)*sqrt(U1_reference);
R0 = I_intval(I_sup(R0_enclosure));

transport_gap = L2_target-U1_reference;
if ~(I_inf(transport_gap) > 0)
    error('calc_ddlami_residual_bounds:TransportGapNotCertified', ...
        ['The affine eigenspace transport requires ', ...
         'L_2(target)>U_1(reference).']);
end

eta_enclosure = L2_target/transport_gap*R0;
eta = I_intval(I_sup(eta_enclosure));
s2_enclosure = eta^2/L2_target;
s2_upper = I_sup(s2_enclosure);
if ~isfinite(s2_upper) || s2_upper < 0 || ~(s2_upper < 1)
    error('calc_ddlami_residual_bounds:TransportAngleNotCertified', ...
        ['The affine eigenspace transport requires ', ...
         'eta^2/L_2(target)<1.']);
end
s2 = I_intval(max(s2_upper,0));

cos_transport = sqrt(I_intval(1)-s2);
cos_transport_lower = I_inf(cos_transport);
phase_loss = I_intval(max(I_sup( ...
    I_intval(1)-I_intval(cos_transport_lower)),0));

dE_sq = eta^2+U1_target*phase_loss^2;
dE = I_intval(I_sup(sqrt( ...
    I_intval(max(I_sup(dE_sq),0)))));
d0_sq = 2*phase_loss;
d0 = I_intval(I_sup(sqrt( ...
    I_intval(max(I_sup(d0_sq),0)))));

eps_a_enclosure = sqrt(qmax)*eps_a_reference+dE;
eps_0_enclosure = eps_0_reference+d0;
eps_a = I_intval(I_sup(eps_a_enclosure));
eps_0 = I_intval(I_sup(eps_0_enclosure));

info = struct();
info.used = true;
info.Q = Q;
info.E = E;
info.E_norm = E_norm;
info.qmin = qmin;
info.qmax = qmax;
info.U1_reference = U1_reference;
info.L2_target = L2_target;
info.U1_target = U1_target;
info.transport_gap = transport_gap;
info.R0 = R0;
info.eta = eta;
info.s2 = s2;
info.cos_transport_lower = cos_transport_lower;
info.dE = dE;
info.d0 = d0;
info.eps_a_reference = eps_a_reference;
info.eps_0_reference = eps_0_reference;
info.eps_a = eps_a;
info.eps_0 = eps_0;
end


function A = local_shape_matrix(P, A_xx, A_xy, A_yy)
A = P(1,1)*A_xx+(P(1,2)+P(2,1))*A_xy+P(2,2)*A_yy;
end


function tf = local_has_interval_width(x)
tf = false;
for k = 1:length(x)
    if I_inf(x(k)) < I_sup(x(k))
        tf = true;
        return;
    end
end
end


function [lower_bound,info] = local_lambda2_analytic_lower( ...
    option,x_triangle,y_triangle,lambda1_bounds)
% Certified analytic lower bound for lambda_2 of
% conv{(0,0),(1,0),(x,y)}.
%
% Krahn--Szego:
%   lambda_2 >= 2*pi*j_{0,1}^2/area
%            >= 4*pi*(2.4048)^2/y.
%
% Lu--Rowlett (fundamental gap of a triangle):
%   lambda_2-lambda_1 >= 64*pi^2/(9*d^2).

if islogical(option)
    if ~option
        lower_bound = I_intval(0);
        info = struct('used',false,'method','disabled');
        return;
    end
    use_standard_bounds = true;
elseif ischar(option)
    normalized = lower(strtrim(option));
    if any(strcmp(normalized,{'none','off','false'}))
        lower_bound = I_intval(0);
        info = struct('used',false,'method','disabled');
        return;
    elseif any(strcmp(normalized, ...
            {'auto','analytic','krahn-szego','lu-rowlett'}))
        use_standard_bounds = true;
    else
        error('calc_ddlami_residual_bounds:BadAnalyticLambda2Option', ...
            'Unrecognized lambda2_analytic_lower option.');
    end
else
    if numel(option) ~= 1
        error('calc_ddlami_residual_bounds:BadAnalyticLambda2Option', ...
            'An explicit analytic lambda_2 lower bound must be scalar.');
    end
    explicit_lower = I_inf(I_intval(option));
    if ~isfinite(explicit_lower) || explicit_lower < 0
        error('calc_ddlami_residual_bounds:BadAnalyticLambda2Option', ...
            ['An explicit analytic lambda_2 lower bound must be ', ...
             'finite and nonnegative.']);
    end
    lower_bound = I_intval(explicit_lower);
    info = struct( ...
        'used',true,'method','user-supplied', ...
        'combined_lower',lower_bound);
    return;
end

if ~use_standard_bounds || ~(I_inf(y_triangle) > 0)
    error('calc_ddlami_residual_bounds:BadAnalyticLambda2Geometry', ...
        'The analytic lambda_2 bounds require a positive triangle height.');
end

j01_lower = I_intval('2.4048');
krahn_enclosure = 4*I_pi*j01_lower^2/y_triangle;
krahn_lower = I_intval(I_inf(krahn_enclosure));

side_left = sqrt(x_triangle^2+y_triangle^2);
side_right = sqrt((I_intval(1)-x_triangle)^2+y_triangle^2);
diameter_upper = I_intval(max([1, ...
    I_sup(side_left),I_sup(side_right)]));
gap_enclosure = 64*I_pi^2/(9*diameter_upper^2);
gap_lower = I_intval(I_inf(gap_enclosure));
lu_rowlett_enclosure = ...
    I_intval(I_inf(lambda1_bounds))+gap_lower;
lu_rowlett_lower = I_intval(I_inf(lu_rowlett_enclosure));

lower_bound = I_intval(max( ...
    I_inf(krahn_lower),I_inf(lu_rowlett_lower)));

info = struct();
info.used = true;
info.method = 'max(Krahn-Szego,Lu-Rowlett)';
info.j01_lower = j01_lower;
info.krahn_lower = krahn_lower;
info.side_left_upper = I_intval(I_sup(side_left));
info.side_right_upper = I_intval(I_sup(side_right));
info.diameter_upper = diameter_upper;
info.fundamental_gap_lower = gap_lower;
info.lambda1_lower = I_intval(I_inf(lambda1_bounds));
info.lu_rowlett_lambda2_lower = lu_rowlett_lower;
info.combined_lower = lower_bound;
end


function strengthened = ...
    local_raise_interval_lower(enclosure,new_lower,label)
global INTERVAL_MODE
lo = max(I_inf(enclosure),I_inf(new_lower));
hi = I_sup(enclosure);
if lo > hi
    if INTERVAL_MODE
        error('calc_ddlami_residual_bounds:InconsistentAnalyticLowerBound', ...
            ['The analytic lower bound for %s exceeds its certified ', ...
             'upper endpoint.'],label);
    else
        % INTERVAL_MODE=0 is an assembly smoke mode: the legacy
        % lower-bound routine may return a non-enclosing point value.
        strengthened = lo;
        return;
    end
end
strengthened = I_infsup(lo,hi);
end


function enclosure = ...
    local_interval_from_lower_upper(lower_bound,upper_bound,label)
lo = I_inf(lower_bound);
hi = I_sup(upper_bound);
if lo > hi
    error('calc_ddlami_residual_bounds:InconsistentEigenvalueBounds', ...
        ['The certified lower bound for %s exceeds its conforming ', ...
         'Rayleigh upper bound.'],label);
end
enclosure = I_infsup(lo,hi);
end


function [inside, boundary, ndof] = local_dirichlet_dofs(mesh, order)
nv = size(mesh.nodes,1);
ne = size(mesh.edges,1);
nt = size(mesh.elements,1);
ndof = nv+(order-1)*ne+(order-1)*(order-2)/2*nt;

boundary_vertices = unique(mesh.boundary_edges(:));
boundary = boundary_vertices(:);
if order > 1
    [found, boundary_edge_ids] = ismember( ...
        sort(mesh.boundary_edges,2),sort(mesh.edges,2),'rows');
    if ~all(found)
        error('calc_ddlami_residual_bounds:BoundaryEdgeMissing', ...
            'A boundary edge is absent from mesh.edges.');
    end
    edge_dofs = [];
    for k = 1:numel(boundary_edge_ids)
        edge_id = boundary_edge_ids(k);
        edge_dofs = [edge_dofs, ...
            nv+(order-1)*(edge_id-1)+(1:(order-1))]; %#ok<AGROW>
    end
    boundary = unique([boundary(:);edge_dofs(:)]);
end
inside = (1:ndof)';
inside(boundary) = [];
end


function [factor_lower,factor_upper] = ...
    local_affine_metric_eigenvalue_factors(S_inv)
% Rigorous extrema for eigenvalues of
%   [1 alpha;0 beta]*[1 alpha;0 beta]'.
% The eigenvalues are even in alpha.  The minimum/maximum over the narrow
% affine cell occur among |alpha|=0 or max|alpha| and the two beta
% endpoints.  Each scalar candidate is evaluated in interval arithmetic so
% that the selected endpoints remain outward rounded.
alphaI = S_inv(1,2);
betaI = S_inv(2,2);
if ~(I_inf(betaI) > 0)
    error('calc_ddlami_residual_bounds:BadAffineMetric', ...
        'The affine height ratio is not certified positive.');
end

alpha_max = max(abs(I_inf(alphaI)),abs(I_sup(alphaI)));
beta_min = min(I_inf(betaI),I_sup(betaI));
beta_max = max(I_inf(betaI),I_sup(betaI));

lower_candidate = Inf;
upper_candidate = -Inf;
for alpha = [0,alpha_max]
    for beta = [beta_min,beta_max]
        [emin,emax] = local_triangular_metric_eigs(alpha,beta);
        lower_candidate = min(lower_candidate,I_inf(emin));
        upper_candidate = max(upper_candidate,I_sup(emax));
    end
end
factor_lower = I_intval(max(lower_candidate,0));
factor_upper = I_intval(upper_candidate);
end


function [emin,emax] = local_triangular_metric_eigs(alpha,beta)
alpha = I_intval(alpha);
beta = I_intval(beta);
trace_metric = I_intval(1)+alpha^2+beta^2;
det_metric = beta^2;
disc_raw = trace_metric^2-4*det_metric;
disc = I_infsup(max(I_inf(disc_raw),0),max(I_sup(disc_raw),0));
root_disc = sqrt(disc);
emin = (trace_metric-root_disc)/2;
emax = (trace_metric+root_disc)/2;
end


function scaled = local_scale_positive_eigen_bounds( ...
    bounds,factor_lower,factor_upper)
scaled = I_intval(zeros(size(bounds)));
for k = 1:length(bounds)
    lo = I_inf(factor_lower*bounds(k));
    hi = I_sup(factor_upper*bounds(k));
    scaled(k) = I_infsup(lo,hi);
end
end


function [Axx_target,Axy_target,Ayy_target] = ...
    local_transport_gradient_components( ...
    Axx_reference,Axy_reference,Ayy_reference,S_inv,detS)
% If grad_target^T=grad_reference^T*S_inv, its x/y components are
% determined by the two columns of S_inv.
qx = S_inv(:,1);
qy = S_inv(:,2);
Axx_target = detS*( ...
    qx(1)^2*Axx_reference ...
    +2*qx(1)*qx(2)*Axy_reference ...
    +qx(2)^2*Ayy_reference);
Ayy_target = detS*( ...
    qy(1)^2*Axx_reference ...
    +2*qy(1)*qy(2)*Axy_reference ...
    +qy(2)^2*Ayy_reference);
Axy_target = detS*( ...
    qx(1)*qy(1)*Axx_reference ...
    +(qx(1)*qy(2)+qx(2)*qy(1))*Axy_reference ...
    +qx(2)*qy(2)*Ayy_reference);
end


function covariance = local_reference_gradient_covariance( ...
    z,w,Axx_reference,Axy_reference,Ayy_reference)
covariance = [z'*Axx_reference*w; ...
              z'*Axy_reference*w; ...
              z'*Ayy_reference*w];
end


function value = local_transported_covariance_form( ...
    covariance,P,S_inv,detS)
% det(S) * grad_ref(z)' * S^{-1} P S^{-T} * grad_ref(w).
metric = S_inv*P*S_inv';
value = detS*( ...
    metric(1,1)*covariance(1) ...
    +2*metric(1,2)*covariance(2) ...
    +metric(2,2)*covariance(3));
end


function value = local_covariance_metric_form(covariance,metric)
% grad_ref(z)' * metric * grad_ref(w), integrated on the reference mesh.
value = metric(1,1)*covariance(1) ...
    +2*metric(1,2)*covariance(2) ...
    +metric(2,2)*covariance(3);
end


function value = local_interval_abs_upper(x)
value = I_intval(max(abs(I_inf(x)),abs(I_sup(x))));
end


function mesh = local_affine_map_triangle_mesh( ...
    mesh_reference,x_reference,y_reference,x_target,y_target)
% Interval-safe counterpart of affine_map_triangle_mesh.  Constructing the
% domain array entry-by-entry also works with INTLAB under Octave, whereas
% mixed double/intval concatenation is MATLAB-specific.
mesh = mesh_reference;
alpha = (x_target-x_reference)/y_reference;
beta = y_target/y_reference;
xy = mesh_reference.nodes;
mapped_x = xy(:,1)+alpha.*xy(:,2);
mapped_y = beta.*xy(:,2);
mesh.nodes = [mapped_x,mapped_y];
mesh.domain = I_zeros(3,2);
mesh.domain(2,1) = I_intval(1);
mesh.domain(3,1) = x_target;
mesh.domain(3,2) = y_target;
end
