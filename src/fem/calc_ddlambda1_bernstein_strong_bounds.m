function [lambda,dlambda,ddlambda_lower,ddlambda_upper,diagnostics] = ...
    calc_ddlambda1_bernstein_strong_bounds( ...
        triangle,e_direction,degree,spectral_options)
%CALC_DDLAMBDA1_BERNSTEIN_STRONG_BOUNDS  Flux-free lambda_i Hessian bound.
%
% The trial space consists of the interior degree-p Bernstein polynomials
% on one reference triangle.  These global polynomials belong to
% H^2 cap H^1_0, so both the eigenfunction and material-equation residuals
% are ordinary L2 strong residuals.  All polynomial moments are evaluated
% with outward-rounded closed formulas; the Hessian step uses no mesh,
% numerical quadrature, H(div) reconstruction, or verified linear solve.
% Its three spectral endpoints are supplied independently by
% Lehmann--Goerisch, conforming Ritz, and CR--Liu enclosures.
%
% TRIANGLE may contain intervals and must be
% [0,0,1,0,x,y].  A Ritz vector and one bordered material trial are chosen
% at the midpoint.  The same fixed coefficient vectors are valid trial
% functions throughout the whole interval cell.

global INTERVAL_MODE
if isempty(INTERVAL_MODE)
    INTERVAL_MODE = 0;
end
if nargin < 3 || isempty(degree)
    degree = 11;
end
if nargin < 4 || isempty(spectral_options)
    spectral_options = struct();
end
spectral_options = local_spectral_options(spectral_options);
spectral_index = spectral_options.index;

triangle = I_intval(triangle);
e_direction = I_intval(e_direction);
if length(triangle) ~= 6 || length(e_direction) ~= 2
    error('calc_ddlambda1_bernstein_strong_bounds:BadInputShape', ...
        'TRIANGLE has length six and E_DIRECTION has length two.');
end
prefix = [0,0,1,0];
for k = 1:4
    if I_inf(triangle(k)) ~= prefix(k) ...
            || I_sup(triangle(k)) ~= prefix(k)
        error('calc_ddlambda1_bernstein_strong_bounds:NoncanonicalTriangle', ...
            'The triangle must be conv{(0,0),(1,0),(x,y)}.');
    end
end
x = triangle(5);
y = triangle(6);
if ~(I_inf(y) > 0)
    error('calc_ddlambda1_bernstein_strong_bounds:BadHeight', ...
        'The triangle height must be strictly positive.');
end

x0 = I_mid(x);
y0 = I_mid(y);
a = e_direction(1);
b = e_direction(2);

poly = triangle_bernstein_bubble_data(degree);
M = poly.M;

[Q,D,H,P_e,P_ee] = local_metrics(x,y,a,b);
[Q0,D0,~,~,~] = local_metrics( ...
    I_intval(x0),I_intval(y0),a,b);
K0 = local_metric_form(poly,Q0);
Ae0 = local_metric_form(poly,D0);
K = local_metric_form(poly,Q);
Ae = local_metric_form(poly,D);
Aee = local_metric_form(poly,H);

% The floating-point eigensolve only chooses trial coefficients.  Every
% quantity used as a bound is recomputed below with interval arithmetic.
Kd = I_mid(K0);
Md = I_mid(M);
[trial_vectors,trial_values] = eig(Kd,Md);
trial_values = real(diag(trial_values));
[trial_values,order] = sort(trial_values,'ascend');
trial_vectors = real(trial_vectors(:,order));
num_trial_values = spectral_index+1;
if numel(trial_values) < num_trial_values ...
        || any(~isfinite(trial_values(1:num_trial_values))) ...
        || any(trial_values(1:num_trial_values) <= 0)
    error('calc_ddlambda1_bernstein_strong_bounds:BadTrialSpectrum', ...
        ['The Bernstein trial eigenproblem did not return the required ', ...
         'positive values.']);
end

raw_u_double = trial_vectors(:,spectral_index);
mass_double = y0*(raw_u_double'*Md*raw_u_double);
if ~(isfinite(mass_double) && mass_double > 0)
    error('calc_ddlambda1_bernstein_strong_bounds:BadTrialMass', ...
        'The selected Bernstein trial vector has non-positive mass.');
end
u_double = raw_u_double/sqrt(mass_double);

% This interval vector encloses the exactly normalized representative of
% the binary floating-point coefficient vector.
raw_u = I_intval(raw_u_double);
mass_raw = I_intval(y0)*(raw_u'*M*raw_u);
u = raw_u/sqrt(mass_raw);
mass_u = I_intval(y0)*(u'*M*u);

K0d = I_mid(K0);
Ae0d = I_mid(Ae0);
mu_anchor = y0*(u_double'*K0d*u_double);
alpha_anchor = y0*(u_double'*Ae0d*u_double);
m0u = y0*Md*u_double;
bordered = [y0*(K0d-mu_anchor*Md),m0u; m0u',0];
rhs = [-y0*Ae0d*u_double+alpha_anchor*m0u;0];
solution = bordered\rhs;
if any(~isfinite(solution))
    error('calc_ddlambda1_bernstein_strong_bounds:BadBorderedSolve', ...
        'The midpoint bordered trial solve returned a non-finite value.');
end
v_double = solution(1:end-1);
v = I_intval(v_double);

% Fixed-reference scalar forms.  The normalized affine pullback contributes
% y0 rather than the variable target determinant; hence no interval
% normalization factor is repeated in these expressions.
K_uu = I_intval(y0)*(u'*K*u);
Ae_uu = I_intval(y0)*(u'*Ae*u);
Aee_uu = I_intval(y0)*(u'*Aee*u);
mass_uv = I_intval(y0)*(u'*M*v);
mass_vv = I_intval(y0)*(v'*M*v);
K_vv = I_intval(y0)*(v'*K*v);
Ae_uv = I_intval(y0)*(u'*Ae*v);

mu = K_uu/mass_u;
alpha = Ae_uu/mass_u;
forms = struct();
forms.mass_u = mass_u;
forms.alpha_h = alpha;
forms.c_h = Aee_uu/mass_u;
forms.ell_h = -Ae_uv+alpha*mass_uv;
forms.b_h_vv = K_vv-mu*mass_vv;
forms.grad_u_sq = K_uu;
forms.grad_v_sq = K_vv;
forms.l2_v_sq = mass_vv;

% Strong eigenfunction residual at the midpoint, followed by a
% cancellation-free cell perturbation bound.
Eu = poly.E*u;
Ev = poly.E*v;
Dxxu = poly.Dxx*u;
Dxyu = poly.Dxy*u;
Dyyu = poly.Dyy*u;
Dxxv = poly.Dxx*v;
Dxyv = poly.Dxy*v;
Dyyv = poly.Dyy*v;

mu_anchor_I = I_intval(mu_anchor);
alpha_anchor_I = I_intval(alpha_anchor);
eig_hessian_anchor = Q0(1,1)*Dxxu ...
    +(Q0(1,2)+Q0(2,1))*Dxyu+Q0(2,2)*Dyyu;
rho_u_anchor = local_mixed_degree_norm( ...
    poly,eig_hessian_anchor,mu_anchor_I*Eu,y0);

norm_u = local_degree_p_norm(poly,Eu,y0);
norm_v = local_degree_p_norm(poly,Ev,y0);
norm_Dxxu = local_degree_p2_norm(poly,Dxxu,y0);
norm_Dxyu = local_degree_p2_norm(poly,Dxyu,y0);
norm_Dyyu = local_degree_p2_norm(poly,Dyyu,y0);
norm_Dxxv = local_degree_p2_norm(poly,Dxxv,y0);
norm_Dxyv = local_degree_p2_norm(poly,Dxyv,y0);
norm_Dyyv = local_degree_p2_norm(poly,Dyyv,y0);

rho_u = rho_u_anchor ...
    +local_abs_upper(Q(1,1)-Q0(1,1))*norm_Dxxu ...
    +(local_abs_upper(Q(1,2)-Q0(1,2)) ...
      +local_abs_upper(Q(2,1)-Q0(2,1)))*norm_Dxyu ...
    +local_abs_upper(Q(2,2)-Q0(2,2))*norm_Dyyu ...
    +local_abs_upper(mu-mu_anchor_I)*norm_u;
rho_u = I_intval(I_sup(rho_u));

% The eigenvalue endpoints are independent of the Hessian trial and are
% computed on the midpoint triangle before affine transport to the cell.
ritz = verified_ritz_enclosures( ...
    I_intval(trial_vectors(:,1:num_trial_values)),K,M,num_trial_values);
[lambda_bounds,spectral_info] = local_spectral_endpoints( ...
    x,y,x0,y0,mu,ritz,spectral_options);
lambda = lambda_bounds(spectral_index);
lambda_i_upper = I_intval(I_sup(lambda));
lambda_ip1_lower = I_intval(I_inf( ...
    lambda_bounds(spectral_index+1)));
mu_lower = I_intval(I_inf(mu));
mu_upper = I_intval(I_sup(mu));
gap_plus = lambda_ip1_lower-mu_upper;
if spectral_index == 1
    gap_minus = I_intval(Inf);
else
    lambda_im1_upper = I_intval(I_sup( ...
        lambda_bounds(spectral_index-1)));
    gap_minus = mu_lower-lambda_im1_upper;
end
if ~(I_inf(gap_plus) > 0) ...
        || (spectral_index > 1 && ~(I_inf(gap_minus) > 0))
    error('calc_ddlambda1_bernstein_strong_bounds:ResidualGapFailed', ...
        ['The verified adjacent eigenvalue bounds do not separate the ', ...
         'selected Bernstein Rayleigh quotient.']);
end
if spectral_index == 1
    energy_factor_minus = I_intval(0);
else
    energy_factor_minus = lambda_im1_upper/gap_minus^2;
end
residual_gap = I_intval(min(I_inf(gap_minus),I_inf(gap_plus)));

q0 = rho_u/residual_gap;
if ~(I_sup(q0) < 1)
    error('calc_ddlambda1_bernstein_strong_bounds:EigenvectorAngleFailed', ...
        'The strong residual does not certify a singleton eigenvector.');
end
energy_factor_plus = lambda_ip1_lower/gap_plus^2;
energy_factor = max(energy_factor_minus,energy_factor_plus);
qa_sq = rho_u^2*energy_factor;
% The rationalized form avoids cancellation when the trial is very
% accurate: 1-sqrt(1-q0^2) = q0^2/(1+sqrt(1-q0^2)).
phase_loss = q0^2/(I_intval(1)+sqrt(I_intval(1)-q0^2));
eps_0 = I_intval(I_sup(sqrt(2*phase_loss)));
eps_a = I_intval(I_sup(sqrt( ...
    qa_sq+lambda_i_upper*phase_loss^2)));

% The material strong residual is another polynomial norm.  The midpoint
% bordered solve need not be verified: any solve error is part of rho_v.
material_hessian_anchor = D0(1,1)*Dxxu ...
    +(D0(1,2)+D0(2,1))*Dxyu+D0(2,2)*Dyyu ...
    +Q0(1,1)*Dxxv ...
    +(Q0(1,2)+Q0(2,1))*Dxyv+Q0(2,2)*Dyyv;
material_mass_anchor = alpha_anchor_I*Eu+mu_anchor_I*Ev;
rho_v_anchor = local_mixed_degree_norm( ...
    poly,material_hessian_anchor,material_mass_anchor,y0);
rho_v = rho_v_anchor ...
    +local_abs_upper(D(1,1)-D0(1,1))*norm_Dxxu ...
    +(local_abs_upper(D(1,2)-D0(1,2)) ...
      +local_abs_upper(D(2,1)-D0(2,1)))*norm_Dxyu ...
    +local_abs_upper(D(2,2)-D0(2,2))*norm_Dyyu ...
    +local_abs_upper(Q(1,1)-Q0(1,1))*norm_Dxxv ...
    +(local_abs_upper(Q(1,2)-Q0(1,2)) ...
      +local_abs_upper(Q(2,1)-Q0(2,1)))*norm_Dxyv ...
    +local_abs_upper(Q(2,2)-Q0(2,2))*norm_Dyyv ...
    +local_abs_upper(alpha-alpha_anchor_I)*norm_u ...
    +local_abs_upper(mu-mu_anchor_I)*norm_v;
rho_v = I_intval(I_sup(rho_v));

core_data = struct();
core_data.i = spectral_index;
core_data.lambda_bounds = lambda_bounds;
core_data.lambda_h = mu;
core_data.P_e = P_e;
core_data.P_ee = P_ee;
core_data.eps_a = eps_a;
core_data.eps_0 = eps_0;
core_data.scalar_forms = forms;
core_data.material_residual_l2 = rho_v;
[ddlambda_lower,ddlambda_upper,core] = ...
    strong_residual_hessian_enclosure(core_data);

dlambda = I_hull(alpha-core.eps_alpha,alpha+core.eps_alpha);

diagnostics = struct();
diagnostics.estimator = ...
    'signed-strong-residual-with-global-Bernstein-bubbles';
diagnostics.spectral_index = spectral_index;
diagnostics.rigorous = logical(INTERVAL_MODE);
diagnostics.enclosure_scope = ...
    local_scope_name(x,y);
diagnostics.is_uniform_over_target_cell = ...
    local_has_width(x) || local_has_width(y);
diagnostics.target_triangle = triangle;
diagnostics.reference_triangle = [0,0,1,0,x0,y0];
diagnostics.reference_choice = 'target-cell-midpoint';
diagnostics.high_mode_truncation = false;
diagnostics.explicit_eigenfunctions = false;
diagnostics.Hdiv_flux_reconstruction = false;
diagnostics.mesh_used = spectral_info.mesh_used;
diagnostics.hessian_mesh_used = false;
diagnostics.spectral_mesh_used = spectral_info.mesh_used;
diagnostics.numerical_quadrature_used = false;
diagnostics.polynomial_degree = degree;
diagnostics.trial_dimension = poly.dimension;
diagnostics.lambda_bounds = lambda_bounds;
diagnostics.lambda_h = mu;
diagnostics.spectral_endpoints = spectral_info;
diagnostics.lambda_i_LG_lower = ...
    spectral_info.lambda_LG_target_lower(spectral_index);
diagnostics.lambda_i_Ritz_upper = ...
    spectral_info.lambda_target_upper(spectral_index);
diagnostics.lambda_ip1_lower = ...
    spectral_info.lambda_next_target_lower;
if spectral_index == 1
    diagnostics.lambda1_LG_lower = ...
        diagnostics.lambda_i_LG_lower;
    diagnostics.lambda1_Ritz_upper = ...
        diagnostics.lambda_i_Ritz_upper;
    diagnostics.lambda2_CR_lower = ...
        diagnostics.lambda_ip1_lower;
end
diagnostics.residual_gap = residual_gap;
diagnostics.gap_minus = gap_minus;
diagnostics.gap_plus = gap_plus;
diagnostics.energy_factor_minus = energy_factor_minus;
diagnostics.energy_factor_plus = energy_factor_plus;
diagnostics.trial_ritz = ritz;
diagnostics.eps_a = eps_a;
diagnostics.eps_0 = eps_0;
diagnostics.eigenfunction_strong_residual = rho_u;
diagnostics.material_strong_residual = rho_v;
diagnostics.eigenfunction_residual_anchor = rho_u_anchor;
diagnostics.material_residual_anchor = rho_v_anchor;
diagnostics.u_coefficients = u;
diagnostics.v_coefficients = v;
diagnostics.P_e = P_e;
diagnostics.P_ee = P_ee;
diagnostics.Q = Q;
diagnostics.D = D;
diagnostics.H = H;
diagnostics.eig_factor_lower = spectral_info.metric_factor_lower;
diagnostics.eig_factor_upper = spectral_info.metric_factor_upper;
diagnostics.flux = struct( ...
    'skipped',true, ...
    'reason','global Bernstein trials make the strong residual an L2 function');
diagnostics.core = core;
diagnostics.lami = lambda;
diagnostics.dlami = dlambda;
diagnostics.ddlam_lower = ddlambda_lower;
diagnostics.ddlam_upper = ddlambda_upper;
end


function [Q,D,H,P_e,P_ee] = local_metrics(x,y,a,b)
oneI = I_intval(1);
zeroI = I_intval(0);
F_inv = [oneI,-x/y;zeroI,oneI/y];
P_e = [zeroI,-a/y;-a/y,-2*b/y];
P_ee = [2*a^2/y^2,4*a*b/y^2; ...
        4*a*b/y^2,6*b^2/y^2];
Q = F_inv*F_inv';
D = F_inv*P_e*F_inv';
H = F_inv*P_ee*F_inv';
end


function A = local_metric_form(poly,metric)
A = metric(1,1)*poly.Axx ...
    +metric(1,2)*poly.Axy ...
    +metric(2,1)*poly.Axy' ...
    +metric(2,2)*poly.Ayy;
A = I_hull(A,A');
end


function value = local_mixed_degree_norm(poly,c2,cp,y0)
sq = c2'*poly.Gp2p2*c2 ...
    +2*(c2'*poly.Gp2p*cp)+cp'*poly.Gpp*cp;
upper = I_sup(I_intval(y0)*sq);
if ~isfinite(upper) || upper < -1e-12
    error('calc_ddlambda1_bernstein_strong_bounds:BadResidualNorm', ...
        'A polynomial residual norm is not finite and nonnegative.');
end
value = I_intval(I_sup(sqrt(I_intval(max(upper,0)))));
end


function value = local_degree_p2_norm(poly,c2,y0)
sq = c2'*poly.Gp2p2*c2;
value = local_sqrt_upper(I_intval(y0)*sq);
end


function value = local_degree_p_norm(poly,cp,y0)
sq = cp'*poly.Gpp*cp;
value = local_sqrt_upper(I_intval(y0)*sq);
end


function value = local_sqrt_upper(x)
upper = I_sup(x);
if ~isfinite(upper) || upper < -1e-12
    error('calc_ddlambda1_bernstein_strong_bounds:BadPolynomialNorm', ...
        'A polynomial norm is not finite and nonnegative.');
end
value = I_intval(I_sup(sqrt(I_intval(max(upper,0)))));
end


function options = local_spectral_options(options)
if ~isstruct(options)
    error('calc_ddlambda1_bernstein_strong_bounds:BadSpectralOptions', ...
        'SPECTRAL_OPTIONS must be a struct.');
end
defaults = struct('N_LG',16,'N_rho',64,'bound_order',2,'index',1);
names = fieldnames(defaults);
for k = 1:numel(names)
    name = names{k};
    if ~isfield(options,name) || isempty(options.(name))
        options.(name) = defaults.(name);
    end
    value = options.(name);
    if ~isscalar(value) || ~isfinite(value) ...
            || value < 1 || value ~= floor(value)
        error('calc_ddlambda1_bernstein_strong_bounds:BadSpectralOptions', ...
            'SPECTRAL_OPTIONS.%s must be a positive integer.',name);
    end
end
end


function [bounds,info] = local_spectral_endpoints( ...
    x,y,x0,y0,mu,bernstein_ritz,options)
global INTERVAL_MODE
index = options.index;
if isfield(options,'reference_certificate') ...
        && ~isempty(options.reference_certificate)
    certificate = options.reference_certificate;
    required = {'schema','reference_triangle','rigorous'};
    for k = 1:numel(required)
        if ~isfield(certificate,required{k})
            error( ...
                'calc_ddlambda1_bernstein_strong_bounds:BadSpectralCertificate', ...
                'Missing reference_certificate.%s.',required{k});
        end
    end
    if ~strcmp(char(certificate.schema), ...
            'lowerboundsineq.triangle-spectral-cr-lg-certificate.v1')
        error( ...
            'calc_ddlambda1_bernstein_strong_bounds:BadSpectralCertificate', ...
            'Unsupported reference-certificate schema.');
    end
    [L_reference,U_reference,Lnext_reference,next_source] = ...
        local_certificate_endpoints(certificate,index);
    if ~local_boolean_scalar(certificate.rigorous)
        error( ...
            'calc_ddlambda1_bernstein_strong_bounds:BadSpectralCertificate', ...
            'reference_certificate.rigorous must be scalar 0 or 1.');
    end
    certificate_rigorous = local_true_scalar(certificate.rigorous);
    if INTERVAL_MODE && ~certificate_rigorous
        error( ...
            'calc_ddlambda1_bernstein_strong_bounds:UnverifiedCertificate', ...
            ['A supplied spectral certificate must set rigorous=true ', ...
             'in interval mode.']);
    end
    spectral_reference = certificate;
    mesh_used = isfield(certificate,'mesh_used') ...
        && logical(certificate.mesh_used);
    method = 'supplied-LG-Ritz-CR-certificate';
    reference_triangle = I_intval(certificate.reference_triangle);
    if length(reference_triangle) ~= 6
        error( ...
            'calc_ddlambda1_bernstein_strong_bounds:BadSpectralCertificate', ...
            'reference_certificate.reference_triangle has length six.');
    end
    prefix = [0,0,1,0];
    for k = 1:4
        if I_inf(reference_triangle(k)) ~= prefix(k) ...
                || I_sup(reference_triangle(k)) ~= prefix(k)
            error( ...
                'calc_ddlambda1_bernstein_strong_bounds:BadSpectralCertificate', ...
                'The spectral reference triangle is not canonical.');
        end
    end
    x_reference = reference_triangle(5);
    y_reference = reference_triangle(6);
    if local_has_width(x_reference) || local_has_width(y_reference) ...
            || ~(I_inf(y_reference) > 0)
        error( ...
            'calc_ddlambda1_bernstein_strong_bounds:BadSpectralCertificate', ...
            ['The spectral reference must be a point triangle of ', ...
             'positive height.']);
    end
else
    x_reference = I_intval(x0);
    y_reference = I_intval(y0);
    reference_triangle = ...
        I_intval([0,0,1,0,x_reference,y_reference]);
    spectral_reference = triangle_spectral_certificate_cr_lg( ...
        reference_triangle,options.N_LG, ...
        options.N_rho,options.bound_order,index);
    [L_reference,U_reference,Lnext_reference,next_source] = ...
        local_certificate_endpoints(spectral_reference,index);
    certificate_rigorous = logical(spectral_reference.rigorous);
    mesh_used = logical(spectral_reference.mesh_used);
    method = spectral_reference.method;
end

if any(I_inf(L_reference) <= 0) ...
        || any(I_sup(U_reference) < I_inf(L_reference)) ...
        || ~(I_inf(Lnext_reference) > I_sup(U_reference(index))) ...
        || (index > 1 ...
            && ~(I_inf(L_reference(index)) ...
                 > I_sup(U_reference(index-1))))
    error('calc_ddlambda1_bernstein_strong_bounds:BadSpectralCertificate', ...
        ['The reference endpoints do not certify the selected ', ...
         'eigenvalue as simple.']);
end

S_inv = [I_intval(1),(x_reference-x)/y; ...
         I_intval(0),y_reference/y];
[qmin,qmax] = local_affine_metric_factors(S_inv);
L_target = I_intval(zeros(index,1));
U_transport = I_intval(zeros(index,1));
U_target = I_intval(zeros(index,1));
for k = 1:index
    L_target(k) = I_intval(I_inf(qmin*L_reference(k)));
    U_transport(k) = I_intval(I_sup(qmax*U_reference(k)));
    U_target(k) = I_intval(min( ...
        I_sup(U_transport(k)),I_sup(bernstein_ritz(k))));
end
if index == 1
    U_target(1) = I_intval(min(I_sup(U_target(1)),I_sup(mu)));
end
Lnext_target = I_intval(I_inf(qmin*Lnext_reference));
Unext_target = I_intval(I_sup(bernstein_ritz(index+1)));
if any(I_inf(L_target) <= 0) ...
        || any(I_inf(L_target) > I_sup(U_target)) ...
        || I_inf(Lnext_target) > I_sup(Unext_target) ...
        || ~(I_inf(Lnext_target) > I_sup(U_target(index))) ...
        || (index > 1 ...
            && ~(I_inf(L_target(index)) ...
                 > I_sup(U_target(index-1))))
    error('calc_ddlambda1_bernstein_strong_bounds:TransportedSpectrumInvalid', ...
        ['The transported endpoints do not certify the selected ', ...
         'eigenvalue as simple.']);
end

bounds = I_intval(zeros(index+1,1));
for k = 1:index
    bounds(k) = I_infsup(I_inf(L_target(k)),I_sup(U_target(k)));
end
bounds(index+1) = ...
    I_infsup(I_inf(Lnext_target),I_sup(Unext_target));
info = struct();
info.method = method;
info.rigorous = certificate_rigorous;
info.index = index;
info.mesh_used = mesh_used;
info.reference_triangle = reference_triangle;
info.reference_certificate = spectral_reference;
info.S_inv = S_inv;
info.metric_factor_lower = qmin;
info.metric_factor_upper = qmax;
info.lambda_LG_reference_lower = L_reference;
info.lambda_Ritz_reference_upper = U_reference;
info.lambda_next_reference_lower = Lnext_reference;
info.lambda_next_reference_source = next_source;
info.lambda_LG_target_lower = L_target;
info.lambda_transport_upper = U_transport;
info.lambda_target_upper = U_target;
info.lambda_next_target_lower = Lnext_target;
info.lambda_next_Bernstein_Ritz_upper = Unext_target;
if index == 1
    info.lambda1_LG_reference_lower = L_reference(1);
    info.lambda1_Ritz_reference_upper = U_reference(1);
    info.lambda2_CR_reference_lower = Lnext_reference;
    info.lambda1_LG_target_lower = L_target(1);
    info.lambda1_transport_upper = U_transport(1);
    info.lambda1_target_upper = U_target(1);
    info.lambda2_CR_target_lower = Lnext_target;
    info.lambda2_Bernstein_Ritz_upper = Unext_target;
end
end


function [L,U,Lnext,next_source] = ...
    local_certificate_endpoints(certificate,index)
if isfield(certificate,'lambda_LG_lower') ...
        && isfield(certificate,'lambda_Ritz_upper')
    L_all = I_intval(certificate.lambda_LG_lower(:));
    U_all = I_intval(certificate.lambda_Ritz_upper(:));
elseif index == 1 ...
        && isfield(certificate,'lambda1_LG_lower') ...
        && isfield(certificate,'lambda1_Ritz_upper')
    L_all = I_intval(certificate.lambda1_LG_lower);
    U_all = I_intval(certificate.lambda1_Ritz_upper);
else
    error('calc_ddlambda1_bernstein_strong_bounds:BadSpectralCertificate', ...
        'The certificate lacks indexed LG/Ritz endpoints.');
end
if length(L_all) < index || length(U_all) < index
    error('calc_ddlambda1_bernstein_strong_bounds:BadSpectralCertificate', ...
        'The certificate does not reach the requested spectral index.');
end
L = I_intval(I_inf(L_all(1:index)));
U = I_intval(I_sup(U_all(1:index)));
if length(L_all) >= index+1
    Lnext = I_intval(I_inf(L_all(index+1)));
    next_source = 'LG-lower-endpoint';
elseif isfield(certificate,'lambda_next_CR_lower')
    Lnext = I_intval(I_inf(certificate.lambda_next_CR_lower));
    next_source = 'CR-Liu-lower-endpoint';
elseif index == 1 && isfield(certificate,'lambda2_CR_lower')
    Lnext = I_intval(I_inf(certificate.lambda2_CR_lower));
    next_source = 'CR-Liu-lower-endpoint';
else
    error('calc_ddlambda1_bernstein_strong_bounds:BadSpectralCertificate', ...
        'The certificate lacks a lower endpoint for lambda_{i+1}.');
end
end


function [factor_lower,factor_upper] = ...
    local_affine_metric_factors(S_inv)
alphaI = S_inv(1,2);
betaI = S_inv(2,2);
if ~(I_inf(betaI) > 0)
    error('calc_ddlambda1_bernstein_strong_bounds:BadAffineMetric', ...
        'The relative affine height is not certified positive.');
end
alpha_max = max(abs(I_inf(alphaI)),abs(I_sup(alphaI)));
beta_min = I_inf(betaI);
beta_max = I_sup(betaI);
lower = Inf;
upper = -Inf;
for alpha = [0,alpha_max]
    for beta = [beta_min,beta_max]
        traceQ = I_intval(1)+I_intval(alpha)^2+I_intval(beta)^2;
        detQ = I_intval(beta)^2;
        disc0 = traceQ^2-4*detQ;
        disc = I_infsup(max(I_inf(disc0),0),max(I_sup(disc0),0));
        root = sqrt(disc);
        lower = min(lower,I_inf((traceQ-root)/2));
        upper = max(upper,I_sup((traceQ+root)/2));
    end
end
factor_lower = I_intval(max(lower,0));
factor_upper = I_intval(upper);
end


function value = local_abs_upper(x)
value = I_intval(max(abs(I_inf(x)),abs(I_sup(x))));
end


function tf = local_has_width(x)
tf = I_inf(x) ~= I_sup(x);
end


function tf = local_true_scalar(value)
tf = local_boolean_scalar(value) ...
    && double(value) == 1;
end


function tf = local_boolean_scalar(value)
tf = isscalar(value) ...
    && (islogical(value) || isnumeric(value)) ...
    && isfinite(double(value)) ...
    && any(double(value) == [0,1]);
end


function name = local_scope_name(x,y)
if local_has_width(x) || local_has_width(y)
    name = 'uniform-over-target-affine-cell';
else
    name = 'pointwise-target-triangle';
end
end
