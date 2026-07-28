function [verified,J_lower,diagnostics] = ...
    verify_J_both_positive(cell_data)
%VERIFY_J_BOTH_POSITIVE  One spectral solve for both Omega_mid functionals.
%
% The eigenvalue lower bound and geometric interval data are common to J1
% and J2.  A CR-first cell escalates to the strict LG implementation if
% either functional is not certified.  An LG-marked cell starts directly
% with LG.  Every accepted lower bound must be finite and positive.

[area_bounds,perimeter_bounds] = compute_geometry_bounds( ...
    cell_data.x_inf,cell_data.x_sup, ...
    cell_data.theta_inf,cell_data.theta_sup);

if logical(local_num(cell_data.isLG))
    methods_to_try = 1;
elseif isfield(cell_data,'allow_lg_escalation') ...
        && ~logical(local_num(cell_data.allow_lg_escalation))
    methods_to_try = 0;
else
    methods_to_try = [0,1];
end

J_lower = struct('J1',-Inf,'J2',-Inf);
J_diagnostics = struct('J1',struct(),'J2',struct());
eigenvalue_diagnostics = struct();
lambda1_lower = NaN;
method_used = 'none';
escalated = false;

for attempt = 1:numel(methods_to_try)
    use_lg = methods_to_try(attempt);
    data = cell_data;
    data.neig = 1;
    data.isLG = use_lg;
    [lambda1_lower,eigenvalue_diagnostics] = ...
        cell_lower_eig_bound(data);
    [J_lower.J1,J_diagnostics.J1] = compute_J_lower_bound( ...
        'J1',lambda1_lower(1),area_bounds,perimeter_bounds);
    [J_lower.J2,J_diagnostics.J2] = compute_J_lower_bound( ...
        'J2',lambda1_lower(1),area_bounds,perimeter_bounds);
    method_used = local_ternary(use_lg,'LG','CR');
    escalated = attempt > 1;
    if local_certified(J_lower.J1) && local_certified(J_lower.J2)
        break;
    end
end

verified_J1 = local_certified(J_lower.J1);
verified_J2 = local_certified(J_lower.J2);
verified = verified_J1 && verified_J2;

diagnostics = struct();
diagnostics.verified = verified;
diagnostics.verified_J1 = verified_J1;
diagnostics.verified_J2 = verified_J2;
diagnostics.method_used = method_used;
diagnostics.escalated = escalated;
diagnostics.lambda1_lower = lambda1_lower(1);
diagnostics.area_bounds = area_bounds;
diagnostics.perimeter_bounds = perimeter_bounds;
diagnostics.eigenvalue = eigenvalue_diagnostics;
diagnostics.functionals = J_diagnostics;
end


function tf = local_certified(value)
tf = isfinite(value) && value > 0;
end


function value = local_ternary(condition,yes_value,no_value)
if condition
    value = yes_value;
else
    value = no_value;
end
end


function value = local_num(input)
if ischar(input) || isstring(input)
    value = str2double(input);
else
    value = double(input);
end
end
