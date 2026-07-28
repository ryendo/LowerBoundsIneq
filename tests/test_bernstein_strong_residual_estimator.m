function test_bernstein_strong_residual_estimator(use_current_mode)
%TEST_BERNSTEIN_STRONG_RESIDUAL_ESTIMATOR  Flux-free Hessian regression.

this_file = mfilename('fullpath');
repo_root = fileparts(fileparts(this_file));
addpath(fullfile(repo_root,'src','interval'),'-begin');
addpath(fullfile(repo_root,'src','fem'),'-begin');
addpath(fullfile(repo_root,'src','algorithms'),'-begin');

global INTERVAL_MODE
if nargin < 1 || ~logical(use_current_mode)
    INTERVAL_MODE = 0;
end
rigorous = logical(INTERVAL_MODE);

poly = triangle_bernstein_bubble_data(11);
assert(poly.dimension == 45);
Mmid = I_mid(poly.M);
assert(norm(Mmid-Mmid','fro') < 1e-13);
assert(min(eig(Mmid)) > 0);

% Algebraic one-sided constant: with exact trial eigenfunction data and
% ||g||=2, gap=5, the lower correction is 2*4/5 and the first-eigenvalue
% upper correction is exactly zero.
forms = struct( ...
    'mass_u',1,'alpha_h',0,'c_h',10,'ell_h',0,'b_h_vv',0, ...
    'grad_u_sq',6,'grad_v_sq',0,'l2_v_sq',0);
data = struct( ...
    'i',1, ...
    'lambda_bounds',[5.5;11.0], ...
    'lambda_h',5.5, ...
    'gap_plus_lower',5, ...
    'P_e',zeros(2), ...
    'P_ee',zeros(2), ...
    'eps_a',0, ...
    'eps_0',0, ...
    'scalar_forms',forms, ...
    'material_residual_l2',2);
[lo,hi,info] = strong_residual_hessian_enclosure(data);
assert(abs(I_mid(lo)-8.4) < 1e-10);
assert(abs(I_mid(hi)-10) < 1e-10);
assert(info.upper_residual_correction == 0);
assert(abs(I_mid(info.certified_width_bound)-1.6) < 1e-10);

% Equilateral point: the degree-11 global polynomial certificate should be
% contain the exact first eigenvalue 16*pi^2/3, with all three target
% curvatures carrying their expected signs.  Only the three spectral
% endpoints use FEM: LG for L1, conforming Ritz for U1, and CR--Liu for L2.
y_top = sqrt(3)/2;
spectral_options = struct( ...
    'N_LG',8,'N_rho',24,'bound_order',2);
if ~rigorous
    lambda1_exact = 16*pi^2/3;
    lambda2_exact = 112*pi^2/9;
    spectral_options.reference_certificate = struct( ...
        'schema', ...
            'lowerboundsineq.triangle-spectral-cr-lg-certificate.v1', ...
        'reference_triangle',[0,0,1,0,0.5,y_top], ...
        'lambda1_LG_lower',lambda1_exact*(1-1e-10), ...
        'lambda1_Ritz_upper',lambda1_exact*(1+1e-10), ...
        'lambda2_CR_lower',lambda2_exact*(1-1e-10), ...
        'rigorous',false, ...
        'mesh_used',false);
end
[lambda,dlambda,ddlo,ddhi,diagnostics] = ...
    calc_ddlambda1_bernstein_strong_bounds( ...
        [0,0,1,0,0.5,y_top],[1,0],11,spectral_options);
lambda_exact = 16*pi^2/3;
if rigorous
    assert(I_inf(lambda) <= lambda_exact);
    assert(lambda_exact <= I_sup(lambda));
else
    assert(abs(I_mid(lambda)-lambda_exact) < 1e-6);
end
assert(diagnostics.trial_dimension == 45);
assert(~diagnostics.Hdiv_flux_reconstruction);
assert(~diagnostics.hessian_mesh_used);
assert(diagnostics.spectral_mesh_used == rigorous);
assert(I_inf(diagnostics.lambda2_CR_lower) ...
    > I_sup(diagnostics.lambda1_Ritz_upper));
curv = omega_up_all_functional_curvatures( ...
    'x',0.5,y_top,dlambda,ddlo,ddhi);
assert(curv.J1.bound > 0);
assert(curv.J2.bound > 0);
assert(curv.JUP.bound < 0);

% A production-width cell exercises the uniform metric path.
xI = I_infsup(0.5,0.50025);
yI = I_infsup(y_top-0.00025,y_top);
spectral_options.reference_certificate = ...
    diagnostics.spectral_endpoints.reference_certificate;
[~,dlambda,ddlo,ddhi,diagnostics] = ...
    calc_ddlambda1_bernstein_strong_bounds( ...
        [0,0,1,0,xI,yI],[1,0],11,spectral_options);
curv = omega_up_all_functional_curvatures( ...
    'x',xI,yI,dlambda,ddlo,ddhi);
assert(curv.J1.bound > 0);
assert(curv.J2.bound > 0);
assert(curv.JUP.bound < 0);
if rigorous
    assert(diagnostics.is_uniform_over_target_cell);
else
    assert(strcmp(diagnostics.enclosure_scope,'pointwise-target-triangle'));
end

fprintf('test_bernstein_strong_residual_estimator: PASS\n');
end
