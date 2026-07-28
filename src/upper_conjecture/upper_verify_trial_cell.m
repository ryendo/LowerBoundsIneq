function result = upper_verify_trial_cell(cell_def, mesh, matrices, mode)
%UPPER_VERIFY_TRIAL_CELL Verify one compact cell by a fixed conforming trial.
%
% A conforming Ritz vector is generated at the cell midpoint, frozen, and
% transported affinely from the reference triangle to every T(x,y) in the
% cell.  If E11,E12,E22,M are its four reference integrals, then
%
%   R(x,y)=((y^2+x^2)E11-2xE12+E22)/(y^2 M).
%
% Thus lambda_1(T(x,y)) <= R(x,y) by the min--max principle.  INTLAB
% evaluates this formula and the conjectured upper bound on the full cell.

required = {'id','x_lo','x_hi','y_lo','y_hi'};
for k = 1:numel(required)
    if ~isfield(cell_def, required{k})
        error('Missing cell field "%s".', required{k});
    end
end
if cell_def.x_lo < 0.5 || cell_def.x_hi > 1 || ...
        cell_def.x_lo > cell_def.x_hi || cell_def.y_lo <= 0 || ...
        cell_def.y_lo > cell_def.y_hi
    error('Invalid compact triangle cell %d.', cell_def.id);
end

x_mid = (cell_def.x_lo + cell_def.x_hi)/2;
y_mid = (cell_def.y_lo + cell_def.y_hi)/2;
if isfield(matrices,'backend') && ...
        strcmp(matrices.backend,'lagrange_high_order')
    [coeff,lambda_ritz] = upper_matrix_trial(matrices,x_mid,y_mid);
    invs = upper_matrix_trial_invariants(matrices,coeff,mode);
    backend = sprintf('P%d_N%d',matrices.trial_order,matrices.N_trial);
    mesh_order = matrices.trial_order;
    trial_resolution = matrices.N_trial;
else
    [coeff, lambda_ritz] = upper_midpoint_trial(mesh, matrices, x_mid, y_mid);
    invs = upper_trial_invariants(mesh, coeff, mode);
    backend = sprintf('P1_N%d',mesh.n);
    mesh_order = 1;
    trial_resolution = mesh.n;
end
is_interval = strcmpi(mode, 'interval');

if is_interval
    x = infsup(cell_def.x_lo, cell_def.x_hi);
    y = infsup(cell_def.y_lo, cell_def.y_hi);
    one = intval(1);
    pi_value = intval('pi');
else
    % Double mode is deliberately pointwise and therefore exploratory.
    x = x_mid;
    y = y_mid;
    one = 1;
    pi_value = pi;
end

N = x^2*invs.Exx - 2*x*invs.Exy + invs.Eyy;
rayleigh = invs.Exx/invs.M0 + N/(y^2*invs.M0);
perimeter = one + sqrt(x^2+y^2) + sqrt((one-x)^2+y^2);
target = pi_value^2*perimeter^2/(3*y^2) ...
    + 2*sqrt(3*one)*pi_value^2/(3*y);
natural_margin = target - rayleigh;

% A direct natural interval for target-rayleigh repeats x and y many times
% and becomes needlessly wide near the Omega_up interface.  In certified
% mode use the centered mean-value form instead.  All midpoint quantities
% are point intervals; derivatives are evaluated once on the full cell.
if is_interval
    x0 = intval(x_mid);
    y0 = intval(y_mid);
    N0 = x0^2*invs.Exx - 2*x0*invs.Exy + invs.Eyy;
    R0 = invs.Exx/invs.M0 + N0/(y0^2*invs.M0);
    r0 = sqrt(x0^2+y0^2);
    s0 = sqrt((one-x0)^2+y0^2);
    P0 = one+r0+s0;
    T0 = pi_value^2*P0^2/(3*y0^2) ...
        + 2*sqrt(3*one)*pi_value^2/(3*y0);
    F0 = T0-R0;

    r = sqrt(x^2+y^2);
    s = sqrt((one-x)^2+y^2);
    P = one+r+s;
    Px = x/r + (x-one)/s;
    Py = y/r + y/s;
    Tx = 2*pi_value^2*P*Px/(3*y^2);
    Ty = 2*pi_value^2*P*Py/(3*y^2) ...
        - 2*pi_value^2*P^2/(3*y^3) ...
        - 2*sqrt(3*one)*pi_value^2/(3*y^2);
    Rx = (2*x*invs.Exx-2*invs.Exy)/(y^2*invs.M0);
    Ry = -2*N/(y^3*invs.M0);
    Fx = Tx-Rx;
    Fy = Ty-Ry;
    margin = F0 + Fx*(x-x0) + Fy*(y-y0);
else
    margin = natural_margin;
end

[rayleigh_lo, rayleigh_hi] = upper_interval_endpoints(rayleigh, mode);
[target_lo, target_hi] = upper_interval_endpoints(target, mode);
[margin_lo, margin_hi] = upper_interval_endpoints(margin, mode);
[natural_margin_lo, natural_margin_hi] = ...
    upper_interval_endpoints(natural_margin, mode);

result.cell_id = cell_def.id;
result.mode = lower(mode);
result.rigor = ternary(is_interval, 'certified_interval', 'exploratory_double');
result.verified = is_interval && margin_lo > 0;
result.x_lo = cell_def.x_lo;
result.x_hi = cell_def.x_hi;
result.y_lo = cell_def.y_lo;
result.y_hi = cell_def.y_hi;
result.backend = backend;
result.mesh_order = mesh_order;
result.trial_resolution = trial_resolution;
result.lambda_ritz_midpoint = lambda_ritz;
result.rayleigh_lower = rayleigh_lo;
result.rayleigh_upper = rayleigh_hi;
result.target_lower = target_lo;
result.target_upper = target_hi;
result.margin_lower = margin_lo;
result.margin_upper = margin_hi;
result.natural_margin_lower = natural_margin_lo;
result.natural_margin_upper = natural_margin_hi;
result.enclosure = ternary(is_interval, 'centered_mean_value', 'midpoint_point');
result.coeff_linf = max(abs(coeff));
end

function value = ternary(condition, yes_value, no_value)
if condition
    value = yes_value;
else
    value = no_value;
end
end
