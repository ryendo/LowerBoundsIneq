function [verified, J_lower, diagnostics] = verify_J_positive(conjecture_type, cell_data)
% VERIFY_J_POSITIVE: Verify J >= 0 for a single cell in Omega_mid
%
% Paper Reference: Section 4, Step 2 (Omega_mid verification), Lemma 4.9 / Algorithm 4.
%   Lower bound used (CORRECT vertices):
%       J_k(T^p) >= B_k( p_{i,j} ; lambda_1(T^{p_{i+1,j+1}}) ),
%   geometry at the inner corner p_{i,j}=(x_inf,theta_inf),
%   eigenvalue lower bound at the outer corner p_{i+1,j+1}=(x_sup,theta_sup).
%
% Eigenvalue method (efficiency strategy):
%   "basically CR only; if it fails, also try LG."  The Crouzeix-Raviart lower
%   bound has gap ~ C_h^2 lambda^2, which is small for nearly degenerate triangles
%   (large J there) but exceeds J near the equilateral, where the tight
%   Lehmann-Goerisch (LG) bound is needed. So:
%     * far from the equilateral: try CR first, escalate to LG only if CR fails;
%     * near p0=(1/2,sqrt3/2): start directly with LG (skip the doomed CR attempt).
%
% Inputs:
%   conjecture_type: 'J1' (Laugesen-Siudeja) or 'J2' (Cheeger-type)
%   cell_data: struct with x_inf,x_sup,theta_inf,theta_sup and FEM params
%     mesh_size_lower_cr, (mesh_size_lower_LG, fem_order_lower_LG for LG).
% Outputs:
%   verified   : true iff J_lower > 0
%   J_lower    : rigorous lower bound on J for this cell (best method tried)
%   diagnostics: struct incl. .method_used ('CR' or 'LG') and .escalated

%% Step 1: Geometry bounds (method-independent) -- compute once
x_inf_str = local_str(cell_data.x_inf);
x_sup_str = local_str(cell_data.x_sup);
theta_inf_str = local_str(cell_data.theta_inf);
theta_sup_str = local_str(cell_data.theta_sup);

[area_bounds, perimeter_bounds] = compute_geometry_bounds(x_inf_str, x_sup_str, ...
    theta_inf_str, theta_sup_str);

%% Step 2: CR-first eigenvalue lower bound, escalate to LG if needed
DIST_LG = 0.30;   % cells within this distance of p0 start directly with LG
x_lo = local_num(cell_data.x_inf);  x_hi = local_num(cell_data.x_sup);
th_lo = local_num(cell_data.theta_inf); th_hi = local_num(cell_data.theta_sup);
xc = 0.5*(x_lo + x_hi); thc = 0.5*(th_lo + th_hi); yc = xc*tan(thc);
dist_p0 = hypot(xc - 0.5, yc - sqrt(3)/2);

if dist_p0 < DIST_LG
    methods_to_try = 1;        % near equilateral: LG only
else
    methods_to_try = [0, 1];   % basically CR; escalate to LG on failure
end

J_lower = -inf; method_used = 'none'; escalated = false;
lam1_lower = NaN; J_diag = struct();
for mi = 1:numel(methods_to_try)
    isLG = methods_to_try(mi);
    cd = cell_data; cd.neig = 1; cd.isLG = isLG;

    lam1_lower = cell_lower_eig_bound(cd);
    [J_try, J_diag] = compute_J_lower_bound(conjecture_type, lam1_lower(1), ...
        area_bounds, perimeter_bounds);

    J_lower = J_try;
    if isLG, method_used = 'LG'; else, method_used = 'CR'; end
    if mi > 1, escalated = true; end

    if J_lower > 0
        break;   % verified with the cheapest sufficient method
    end
end

%% Step 3: Verification status
verified = (J_lower > 0);

%% Diagnostics
diagnostics = struct();
diagnostics.conjecture_type = conjecture_type;
diagnostics.cell_data = cell_data;
diagnostics.lam1_lower = lam1_lower(1);
diagnostics.area_bounds = area_bounds;
diagnostics.perimeter_bounds = perimeter_bounds;
diagnostics.J_lower = J_lower;
diagnostics.J_diagnostics = J_diag;
diagnostics.verified = verified;
diagnostics.method_used = method_used;   % 'CR' or 'LG'
diagnostics.escalated = escalated;       % true if CR failed and LG was used
diagnostics.dist_p0 = dist_p0;

end

% ---- helpers ----
function s = local_str(v)
    if ischar(v) || isstring(v), s = char(v); else, s = num2str(v, '%.17g'); end
end

function x = local_num(v)
    if ischar(v) || isstring(v), x = str2double(v); else, x = double(v); end
end
