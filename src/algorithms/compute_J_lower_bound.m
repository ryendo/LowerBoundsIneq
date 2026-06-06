function [J_lower, diagnostics] = compute_J_lower_bound(conjecture_type, lambda1_lower, area_bounds, perimeter_bounds)
% COMPUTE_J_LOWER_BOUND (Lemma-aligned):
% Compute the Lemma lower bound:
%   J_k(T^p) >= B_k( p_{i+1,j} ; lambda1(T^{p_{i+1,j+1}}) )
% They are the two vertex values produced by compute_geometry_bounds:
%   area_bounds(1)      = Area at p_{i+1,j}
%   area_bounds(2)      = Area at p_{i,j+1}
%   perimeter_bounds(1) = Perimeter at p_{i+1,j}
%   perimeter_bounds(2) = Perimeter at p_{i,j+1}
%
% lambda1_lower:
%   - scalar: interpreted as lam_low = lower bound for lambda1(T^{p_{i+1,j+1}})
%   - 2-vector [lam_low, lam_up]: additionally compute Lemma upper bound candidate in diagnostics

    % -----------------------------------------
    % Parse lambda input (keep backward compat)
    % -----------------------------------------
    lam_low = [];
    lam_up  = [];

    if isnumeric(lambda1_lower) && numel(lambda1_lower) == 2
        lam_low = I_intval(lambda1_lower(1));
        lam_up  = I_intval(lambda1_lower(2));
    elseif isstruct(lambda1_lower)
        if isfield(lambda1_lower,'low'), lam_low = I_intval(lambda1_lower.low); end
        if isfield(lambda1_lower,'up'),  lam_up  = I_intval(lambda1_lower.up);  end
        if isempty(lam_low) && isfield(lambda1_lower,'lower'), lam_low = I_intval(lambda1_lower.lower); end
        if isempty(lam_up)  && isfield(lambda1_lower,'upper'), lam_up  = I_intval(lambda1_lower.upper); end
    else
        lam_low = I_intval(lambda1_lower);
    end

    if isempty(lam_low)
        error('lambda1_lower must provide a lower bound (scalar or [low,up] or struct).');
    end

    % -----------------------------------------
    % Lemma vertices geometry (by convention)
    % -----------------------------------------
    % p_{i+1,j} geometry for the LOWER bound
    A_L = I_intval(area_bounds(1));
    P_L = I_intval(perimeter_bounds(1));

    % -----------------------------------------
    % Evaluate B_k at the Lemma lower vertex
    % -----------------------------------------
    if strcmpi(conjecture_type, 'J1')
        % B1(p;Lambda) = Lambda*A - (pi^2/16)*P^2/A - 7*sqrt(3)*pi^2/12
        C1 = I_pi^2 / I_intval('16');
        C2 = I_intval('7') * sqrt(I_intval('3')) * I_pi^2 / I_intval('12');

        term1_L = lam_low * A_L;
        term2_L = C1 * (P_L^2) / A_L;
        J_L = term1_L - term2_L - C2;

        J_lower = I_inf(J_L);

        diagnostics = struct();
        diagnostics.conjecture = 'J1';
        diagnostics.B_eval_point = 'p_{i+1,j}';
        diagnostics.lambda_used = lam_low;
        diagnostics.A_used = A_L;
        diagnostics.P_used = P_L;
        diagnostics.term1 = term1_L;
        diagnostics.term2 = term2_L;
        diagnostics.C2 = C2;
        diagnostics.J_interval = J_L;
        diagnostics.J_lower = J_lower;

    elseif strcmpi(conjecture_type, 'J2')
        % B2(p;Lambda) = Lambda*A - C_* * (P + sqrt(4*pi*A))^2 / (4*A)
        sqrt3 = sqrt(I_intval('3'));
        inner_sqrt = sqrt(I_pi * sqrt3);
        denom = (I_intval('3') + inner_sqrt)^2;
        C_star = I_intval('4') * I_pi^2 / denom;

        term1_L = lam_low * A_L;
        sqrt_term_L = sqrt(I_intval('4') * I_pi * A_L);
        term2_L = C_star * (P_L + sqrt_term_L)^2 / (I_intval('4') * A_L);
        J_L = term1_L - term2_L;

        J_lower = I_inf(J_L);

        diagnostics = struct();
        diagnostics.conjecture = 'J2';
        diagnostics.B_eval_point = 'p_{i+1,j}';
        diagnostics.lambda_used = lam_low;
        diagnostics.A_used = A_L;
        diagnostics.P_used = P_L;
        diagnostics.C_star = C_star;
        diagnostics.term1 = term1_L;
        diagnostics.term2 = term2_L;
        diagnostics.J_interval = J_L;
        diagnostics.J_lower = J_lower;

    else
        error('Unknown conjecture type: %s. Use ''J1'' or ''J2''.', conjecture_type);
    end
end
