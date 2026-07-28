function [ddlam_lower, ddlam_upper, diagnostics] = residual_hessian_enclosure(data)
%RESIDUAL_HESSIAN_ENCLOSURE  Two-sided residual bound for a simple eigenvalue.
%
% This is the algebraic (mesh-independent) part of the shifted-flux
% estimator.  All quantities may be INTLAB intervals.  The required fields
% of DATA are
%
%   i, lambda_bounds, lambda_h, lambda1_lower,
%   u_h, v_h, K, M, A_e, A_ee, P_e, P_ee,
%   eps_a, eps_0, flux_residual.
%
% Here u_h is L2-normalized (an interval normalization is allowed), v_h is
% any conforming approximation to the material derivative, and
%
%   flux_residual >= ||r_h||_{a^*},
%
% where r_h is the full discrete material-equation residual.  An exactly
% equilibrated flux satisfying
%
%   div sigma = -(alpha_h u_h + lambda_h v_h),
%   alpha_h = a_e(u_h,u_h).
%
% gives the special case
% ||r_h||_{a^*} <= ||sigma-(P_e grad u_h+grad v_h)||.  A flux with a
% certified equilibrium defect may instead include its Poincare-weighted
% source residual in flux_residual.
%
% With
%
%   ell_h(z) = -a_e(u_h,z)+alpha_h (u_h,z),
%   b_h(z,z) = a(z,z)-lambda_h (z,z),
%   J_h      = 2 ell_h(v_h)-b_h(v_h,v_h),
%
% the returned enclosure is
%
%   c_h-2J_h-eps_D-2 C_plus R^2 <= lambda_i''
%       <= c_h-2J_h+eps_D+2 C_minus R^2.
%
% The energy-dual resolvent constants are
%
%   C_minus = U_{i-1}/(L_i-U_{i-1})  (zero when i=1),
%   C_plus  = L_{i+1}/(L_{i+1}-U_i).
%
% Consequently the residual correction on the *upper* Hessian bound
% vanishes for the first eigenvalue.

upper_only = isfield(data,'upper_only') && logical(data.upper_only);
global INTERVAL_MODE
if isempty(INTERVAL_MODE)
    INTERVAL_MODE = 0;
end
use_scalar_forms = isfield(data,'scalar_forms');
required = {'i','lambda_bounds','lambda_h','lambda1_lower', ...
    'P_e','P_ee','eps_a','eps_0'};
if use_scalar_forms
    required{end+1} = 'scalar_forms';
else
    required = [required,{'u_h','v_h','K','M','A_e','A_ee'}];
end
if ~upper_only
    required{end+1} = 'flux_residual';
end
for k = 1:numel(required)
    if ~isfield(data, required{k})
        error('residual_hessian_enclosure:MissingField', ...
            'Missing required field data.%s.', required{k});
    end
end

i = data.i;
lams = data.lambda_bounds(:);
if i < 1 || i ~= floor(i)
    error('residual_hessian_enclosure:BadIndex', ...
        'The eigenvalue index i must be a positive integer.');
end
if upper_only && i ~= 1
    error('residual_hessian_enclosure:UpperOnlyNeedsFirstEigenvalue', ...
        'The flux-free upper-only path is valid only for i=1.');
end
if length(lams) < i + 1
    error('residual_hessian_enclosure:NeedNextEigenvalue', ...
        'A certified enclosure of lambda_{i+1} is required.');
end

L_i = I_inf(lams(i));
U_i = I_sup(lams(i));
L_ip1 = I_inf(lams(i+1));
if ~(L_ip1 > U_i)
    error('residual_hessian_enclosure:UpperGapNotCertified', ...
        'Cannot certify lambda_i < lambda_{i+1}.');
end

if i == 1
    C_minus = I_intval(0);
else
    U_im1 = I_sup(lams(i-1));
    if ~(L_i > U_im1)
        error('residual_hessian_enclosure:LowerGapNotCertified', ...
            'Cannot certify lambda_{i-1} < lambda_i.');
    end
    C_minus = I_intval(U_im1) ...
        / (I_intval(L_i)-I_intval(U_im1));
end
C_plus = I_intval(L_ip1)/(I_intval(L_ip1)-I_intval(U_i));

lambda_h = data.lambda_h;
if use_scalar_forms
    forms = data.scalar_forms;
    form_names = {'mass_u','alpha_h','c_h','ell_h','b_h_vv', ...
        'grad_u_sq','grad_v_sq','l2_v_sq'};
    for k = 1:length(form_names)
        if ~isfield(forms,form_names{k})
            error('residual_hessian_enclosure:MissingScalarForm', ...
                'Missing data.scalar_forms.%s.',form_names{k});
        end
    end
    mass_u = forms.mass_u;
    alpha_h = forms.alpha_h;
    c_h = forms.c_h;
    ell_h = forms.ell_h;
    b_h_vv = forms.b_h_vv;
    grad_u = sqrt_nonnegative_upper(forms.grad_u_sq);
    grad_v = sqrt_nonnegative_upper(forms.grad_v_sq);
    l2_v = sqrt_nonnegative_upper(forms.l2_v_sq);
else
    u_h = data.u_h(:);
    v_h = data.v_h(:);
    K = data.K;
    M = data.M;
    A_e = data.A_e;
    A_ee = data.A_ee;
    mass_u = u_h'*M*u_h;
    alpha_h = (u_h'*A_e*u_h)/mass_u;
    c_h = (u_h'*A_ee*u_h)/mass_u;
    ell_h = -u_h'*A_e*v_h+alpha_h*(u_h'*M*v_h);
    b_h_vv = v_h'*K*v_h-lambda_h*(v_h'*M*v_h);
    grad_u = sqrt_nonnegative_upper(u_h'*K*u_h);
    grad_v = sqrt_nonnegative_upper(v_h'*K*v_h);
    l2_v = sqrt_nonnegative_upper(v_h'*M*v_h);
end
if ~(I_inf(mass_u) > 0)
    error('residual_hessian_enclosure:BadNormalization', ...
        'The discrete eigenfunction has non-positive L2 norm.');
end
if INTERVAL_MODE
    normalized = I_inf(mass_u) <= 1 && 1 <= I_sup(mass_u);
else
    normalized = abs(I_mid(mass_u)-1) <= 1e-8;
end
if ~normalized
    error('residual_hessian_enclosure:BadNormalization', ...
        ['u_h must be L2-normalized, and eps_a/eps_0 must refer to ', ...
         'that normalized representative.']);
end

J_h = 2*ell_h-b_h_vv;
center_h = c_h-2*J_h;

p = interval_matrix_norm_upper(data.P_e);
q = interval_matrix_norm_upper(data.P_ee);
eps_a = nonnegative_upper(data.eps_a, 'eps_a');
eps_0 = nonnegative_upper(data.eps_0, 'eps_0');

lambda_difference = lams(i) - lambda_h;
eps_lambda = max(abs(I_inf(lambda_difference)), ...
                 abs(I_sup(lambda_difference)));
eps_lambda = I_intval(eps_lambda);

common_u = 2 * grad_u * eps_a + eps_a^2;
eps_alpha = p * common_u;
eps_c = q * common_u;
abs_alpha_h = interval_abs_upper(alpha_h);

eps_ell = p * eps_a * grad_v ...
    + (eps_alpha + abs_alpha_h * eps_0) * l2_v;
eps_D = eps_c + 4 * eps_ell + 2 * eps_lambda * l2_v^2;

% The shifted flux equilibrates the fully discrete source.  The remaining
% source mismatch is due only to u, alpha and lambda errors.
eps_f = eps_alpha + abs_alpha_h * eps_0 + eps_lambda * l2_v;

lambda1_lower = I_intval(data.lambda1_lower);
if ~(I_inf(lambda1_lower) > 0)
    error('residual_hessian_enclosure:BadPoincareConstant', ...
        'lambda1_lower must be strictly positive.');
end

if upper_only
    eta_flux = I_intval(0);
    R_gradient = I_intval(0);
    R_source = I_intval(0);
    R = I_intval(0);
    beta_opt = I_intval(NaN);
else
    eta_flux = nonnegative_upper(data.flux_residual,'flux_residual');
    R_gradient = eta_flux+p*eps_a;
    R_source = eps_f/sqrt(I_intval(I_inf(lambda1_lower)));

    % This is the beta-optimized form of
    %   (1+beta) R_gradient^2
    %     +(1+1/beta) lambda_1^{-1} eps_f^2.
    R = R_gradient+R_source;
    if I_sup(R_gradient) > 0
        beta_opt = R_source/R_gradient;
    elseif I_sup(R_source) > 0
        beta_opt = I_intval(Inf);
    else
        beta_opt = I_intval(1);
    end
end

if upper_only
    lower_interval = I_intval(-Inf);
else
    lower_interval = center_h-eps_D-2*C_plus*R^2;
end
upper_interval = center_h + eps_D + 2 * C_minus * R^2;

ddlam_lower = I_intval(I_inf(lower_interval));
ddlam_upper = I_intval(I_sup(upper_interval));
if I_inf(ddlam_lower) > I_sup(ddlam_upper)
    error('residual_hessian_enclosure:InvertedBounds', ...
        'Computed lower bound exceeds the upper bound.');
end

diagnostics = struct();
diagnostics.upper_only = upper_only;
diagnostics.used_scalar_forms = use_scalar_forms;
diagnostics.C_minus = C_minus;
diagnostics.C_plus = C_plus;
diagnostics.alpha_h = alpha_h;
diagnostics.c_h = c_h;
diagnostics.ell_h = ell_h;
diagnostics.b_h_vv = b_h_vv;
diagnostics.J_h = J_h;
diagnostics.center_h = center_h;
diagnostics.grad_u_h = grad_u;
diagnostics.grad_v_h = grad_v;
diagnostics.l2_v_h = l2_v;
diagnostics.eps_a = eps_a;
diagnostics.eps_0 = eps_0;
diagnostics.eps_lambda = eps_lambda;
diagnostics.eps_alpha = eps_alpha;
diagnostics.eps_c = eps_c;
diagnostics.eps_ell = eps_ell;
diagnostics.eps_D = eps_D;
diagnostics.eps_f = eps_f;
diagnostics.flux_residual = eta_flux;
diagnostics.R_gradient = R_gradient;
diagnostics.R_source = R_source;
diagnostics.R = R;
diagnostics.R_squared = R^2;
diagnostics.beta_opt = beta_opt;
diagnostics.lower_interval = lower_interval;
diagnostics.upper_interval = upper_interval;
diagnostics.ddlam_lower = ddlam_lower;
diagnostics.ddlam_upper = ddlam_upper;
diagnostics.upper_residual_correction = 2 * C_minus * R^2;
diagnostics.lower_residual_correction = 2 * C_plus * R^2;
end


function value = nonnegative_upper(x, label)
upper = I_sup(x);
if isnan(upper) || isinf(upper) || upper < 0
    error('residual_hessian_enclosure:BadNonnegativeInput', ...
        '%s must have a finite non-negative upper endpoint.', label);
end
value = I_intval(max(upper, 0));
end


function value = sqrt_nonnegative_upper(x)
upper = I_sup(x);
if isnan(upper) || isinf(upper)
    error('residual_hessian_enclosure:NaNNorm', ...
        'A discrete squared norm has a NaN upper endpoint.');
end
root_enclosure = sqrt(I_intval(max(upper,0)));
value = I_intval(I_sup(root_enclosure));
end


function value = interval_abs_upper(x)
value = I_intval(max(abs(I_inf(x)), abs(I_sup(x))));
end


function value = interval_matrix_norm_upper(A)
% INTLAB's norm is interval-valued; retain only its guaranteed upper end.
value = I_intval(I_sup(norm(A, 2)));
end
