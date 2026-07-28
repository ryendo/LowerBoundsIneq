function [ddlam_lower,ddlam_upper,diagnostics] = ...
    strong_residual_hessian_enclosure(data)
%STRONG_RESIDUAL_HESSIAN_ENCLOSURE  Signed L2-residual Hessian bounds.
%
% This is the flux-free algebraic core used with globally H2-conforming
% trial functions.  Let u_h approximate the normalized eigenfunction and
% let v_h be any H2 cap H01 trial function for its material derivative.
% For
%
%   g_h = div(P_e grad u_h + grad v_h)
%         + alpha_h u_h + lambda_h v_h,
%
% DATA.material_residual_l2 must enclose ||g_h||_0.  No H(div) flux
% reconstruction is required.  The signed resolvent gives
%
%  lambda_i'' >= c_h-2J_h-eps_D
%      -2 ( (||g_h||+delta_0)/sqrt(gap_+)
%           +sqrt(C_+) delta_a )^2,
%
%  lambda_i'' <= c_h-2J_h+eps_D
%      +2 ( (||g_h||+delta_0)/sqrt(gap_-)
%           +sqrt(C_-) delta_a )^2,
%
% where the second correction is identically zero for i=1.  The
% eigenpair-induced residual mismatch is split into the gradient part
% delta_a=||P_e|| eps_a and the L2 source mismatch delta_0.  This
% norm-matched split is sharper than first converting delta_0 by Poincare.
%
% Relative to first converting ||g_h||_0 by Poincare and then applying an
% energy-dual resolvent constant, the coefficient of the strong residual
% is reduced from sqrt(C_\pm/lambda_1) to 1/sqrt(gap_\pm).

global INTERVAL_MODE
if isempty(INTERVAL_MODE)
    INTERVAL_MODE = 0;
end

required = {'i','lambda_bounds','lambda_h', ...
    'P_e','P_ee','eps_a','eps_0','scalar_forms', ...
    'material_residual_l2'};
for k = 1:numel(required)
    if ~isfield(data,required{k})
        error('strong_residual_hessian_enclosure:MissingField', ...
            'Missing required field data.%s.',required{k});
    end
end

i = data.i;
if i < 1 || i ~= floor(i)
    error('strong_residual_hessian_enclosure:BadIndex', ...
        'The eigenvalue index must be a positive integer.');
end
lams = I_intval(data.lambda_bounds(:));
if length(lams) < i+1
    error('strong_residual_hessian_enclosure:NeedNextEigenvalue', ...
        'A certified enclosure of lambda_{i+1} is required.');
end

L_i = I_intval(I_inf(lams(i)));
U_i = I_intval(I_sup(lams(i)));
L_ip1 = I_intval(I_inf(lams(i+1)));
if isfield(data,'gap_plus_lower') && ~isempty(data.gap_plus_lower)
    supplied_gap_plus = I_intval(data.gap_plus_lower);
    gap_plus = I_intval(I_inf(supplied_gap_plus));
    gap_plus_source = 'direct-gap-lower';
else
    gap_plus = L_ip1-U_i;
    gap_plus_source = 'adjacent-endpoints';
end
if ~(I_inf(gap_plus) > 0)
    error('strong_residual_hessian_enclosure:UpperGapNotCertified', ...
        'Cannot certify lambda_i < lambda_{i+1}.');
end
if strcmp(gap_plus_source,'direct-gap-lower')
    % lambda_{i+1}/(lambda_{i+1}-lambda_i)
    % = 1+lambda_i/(lambda_{i+1}-lambda_i).
    C_plus = I_intval(1)+U_i/gap_plus;
else
    C_plus = L_ip1/gap_plus;
end

if i == 1
    gap_minus = I_intval(Inf);
    C_minus = I_intval(0);
else
    U_im1 = I_intval(I_sup(lams(i-1)));
    if isfield(data,'gap_minus_lower') && ~isempty(data.gap_minus_lower)
        supplied_gap_minus = I_intval(data.gap_minus_lower);
        gap_minus = I_intval(I_inf(supplied_gap_minus));
        gap_minus_source = 'direct-gap-lower';
    else
        gap_minus = L_i-U_im1;
        gap_minus_source = 'adjacent-endpoints';
    end
    if ~(I_inf(gap_minus) > 0)
        error('strong_residual_hessian_enclosure:LowerGapNotCertified', ...
            'Cannot certify lambda_{i-1} < lambda_i.');
    end
    C_minus = U_im1/gap_minus;
end

forms = data.scalar_forms;
form_names = {'mass_u','alpha_h','c_h','ell_h','b_h_vv', ...
    'grad_u_sq','grad_v_sq','l2_v_sq'};
for k = 1:numel(form_names)
    if ~isfield(forms,form_names{k})
        error('strong_residual_hessian_enclosure:MissingScalarForm', ...
            'Missing data.scalar_forms.%s.',form_names{k});
    end
end

mass_u = I_intval(forms.mass_u);
if ~(I_inf(mass_u) > 0)
    error('strong_residual_hessian_enclosure:BadNormalization', ...
        'The trial eigenfunction has non-positive mass.');
end
if INTERVAL_MODE
    normalized = I_inf(mass_u) <= 1 && 1 <= I_sup(mass_u);
else
    normalized = abs(I_mid(mass_u)-1) <= 1e-8;
end
if ~normalized
    error('strong_residual_hessian_enclosure:BadNormalization', ...
        'The trial eigenfunction must be L2-normalized.');
end

alpha_h = I_intval(forms.alpha_h);
c_h = I_intval(forms.c_h);
ell_h = I_intval(forms.ell_h);
b_h_vv = I_intval(forms.b_h_vv);
J_h = 2*ell_h-b_h_vv;
center_h = c_h-2*J_h;

grad_u = local_sqrt_nonnegative_upper(forms.grad_u_sq);
grad_v = local_sqrt_nonnegative_upper(forms.grad_v_sq);
l2_v = local_sqrt_nonnegative_upper(forms.l2_v_sq);
p = local_matrix_norm_upper(I_intval(data.P_e));
q = local_matrix_norm_upper(I_intval(data.P_ee));
eps_a = local_nonnegative_upper(data.eps_a,'eps_a');
eps_0 = local_nonnegative_upper(data.eps_0,'eps_0');
g_h = local_nonnegative_upper( ...
    data.material_residual_l2,'material_residual_l2');

lambda_h = I_intval(data.lambda_h);
lambda_difference = lams(i)-lambda_h;
eps_lambda = I_intval(max(abs(I_inf(lambda_difference)), ...
                          abs(I_sup(lambda_difference))));

common_u = 2*grad_u*eps_a+eps_a^2;
eps_alpha = p*common_u;
eps_c = q*common_u;
abs_alpha_h = local_abs_upper(alpha_h);
eps_ell = p*eps_a*grad_v ...
    +(eps_alpha+abs_alpha_h*eps_0)*l2_v;
eps_D = eps_c+4*eps_ell+2*eps_lambda*l2_v^2;

eps_f = eps_alpha+abs_alpha_h*eps_0+eps_lambda*l2_v;
delta_a = p*eps_a;
delta_0 = eps_f;

positive_radius = (g_h+delta_0)/sqrt(gap_plus) ...
    +sqrt(C_plus)*delta_a;
lower_correction = 2*positive_radius^2;
if i == 1
    negative_radius = I_intval(0);
    upper_correction = I_intval(0);
else
    negative_radius = (g_h+delta_0)/sqrt(gap_minus) ...
        +sqrt(C_minus)*delta_a;
    upper_correction = 2*negative_radius^2;
end

lower_interval = center_h-eps_D-lower_correction;
upper_interval = center_h+eps_D+upper_correction;
ddlam_lower = I_intval(I_inf(lower_interval));
ddlam_upper = I_intval(I_sup(upper_interval));
if I_inf(ddlam_lower) > I_sup(ddlam_upper)
    error('strong_residual_hessian_enclosure:InvertedBounds', ...
        'Computed lower bound exceeds the upper bound.');
end

diagnostics = struct();
diagnostics.estimator = 'signed-strong-residual-resolvent';
diagnostics.C_minus = C_minus;
diagnostics.C_plus = C_plus;
diagnostics.gap_minus = gap_minus;
diagnostics.gap_plus = gap_plus;
if i == 1
    diagnostics.gap_minus_source = 'empty-below-first-eigenvalue';
else
    diagnostics.gap_minus_source = gap_minus_source;
end
diagnostics.gap_plus_source = gap_plus_source;
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
diagnostics.delta_a = delta_a;
diagnostics.delta_0 = delta_0;
diagnostics.material_residual_l2 = g_h;
diagnostics.positive_radius = positive_radius;
diagnostics.negative_radius = negative_radius;
diagnostics.lower_residual_correction = lower_correction;
diagnostics.upper_residual_correction = upper_correction;
diagnostics.lower_interval = lower_interval;
diagnostics.upper_interval = upper_interval;
diagnostics.ddlam_lower = ddlam_lower;
diagnostics.ddlam_upper = ddlam_upper;
diagnostics.certified_width_bound = ...
    2*eps_D+lower_correction+upper_correction;
end


function value = local_nonnegative_upper(x,name)
x = I_intval(x);
upper = I_sup(x);
if ~isfinite(upper) || I_inf(x) < 0
    error('strong_residual_hessian_enclosure:BadNonnegativeInput', ...
        '%s must be finite and nonnegative.',name);
end
value = I_intval(upper);
end


function value = local_sqrt_nonnegative_upper(x)
x = I_intval(x);
upper = I_sup(x);
if ~isfinite(upper) || upper < 0
    error('strong_residual_hessian_enclosure:BadSquaredNorm', ...
        'A squared norm must have a finite nonnegative upper endpoint.');
end
value = I_intval(I_sup(sqrt(I_intval(max(upper,0)))));
end


function value = local_abs_upper(x)
x = I_intval(x);
value = I_intval(max(abs(I_inf(x)),abs(I_sup(x))));
end


function value = local_matrix_norm_upper(A)
A = I_intval(A);
value = I_intval(I_sup(norm(A,2)));
end
