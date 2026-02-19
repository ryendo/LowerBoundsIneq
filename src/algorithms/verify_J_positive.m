function [verified, J_lower, diagnostics] = verify_J_positive(conjecture_type, cell_data)
% VERIFY_J_POSITIVE: Verify J >= 0 for a single cell in Omega_mid
%
% Paper Reference: Section 4, Step 2 (Omega_mid verification)
%   For each cell in the grid, compute a lower bound for J and verify J >= 0.
%
%
% Inputs:
%   conjecture_type: 'J1' (Laugesen-Siudeja) or 'J2' (Cheeger-type)
%   cell_data: structure with fields:
%     - x_inf, x_sup: x-coordinate bounds (strings or numbers)
%     - theta_inf, theta_sup: theta angle bounds (strings or numbers)
%     - mesh_size_upper, fem_order_upper: FEM parameters for upper bound
%     - mesh_size_lower_cr: mesh size for CR lower bound
%     - isLG: use Lehmann-Goerisch (1) or CR (0)
%     - mesh_size_lower_LG, fem_order_lower_LG: LG parameters (if isLG=1)
%
% Outputs:
%   verified: true if J_lower > 0 (J is positive over the cell)
%   J_lower: rigorous lower bound on J for this cell
%   diagnostics: detailed computation information
%
% Author: Based on paper by R. Endo, X. Liu, P. Mariano
% Date: 2025-01-14

% Add paths
addpath('../FEM_Functions');
addpath('../mode_swith_interface');

%% Step 1: Compute lambda_1 lower bound
% Set neig = 1 for efficiency (only lambda_1 needed for J evaluation)
cell_data_lam1 = cell_data;
cell_data_lam1.neig = 1;

fprintf('  Computing lambda_1 lower bound...\n');
lam1_lower = cell_lower_eig_bound(cell_data_lam1);
fprintf('    lambda_1 lower bound: %.17f\n', lam1_lower(1));

%% Step 2: Compute geometry bounds (area and perimeter)
fprintf('  Computing geometry bounds...\n');

% Convert cell_data to strings if necessary
if isnumeric(cell_data.x_inf)
    x_inf_str = num2str(cell_data.x_inf, '%.17g');
else
    x_inf_str = cell_data.x_inf;
end
if isnumeric(cell_data.x_sup)
    x_sup_str = num2str(cell_data.x_sup, '%.17g');
else
    x_sup_str = cell_data.x_sup;
end
if isnumeric(cell_data.theta_inf)
    theta_inf_str = num2str(cell_data.theta_inf, '%.17g');
else
    theta_inf_str = cell_data.theta_inf;
end
if isnumeric(cell_data.theta_sup)
    theta_sup_str = num2str(cell_data.theta_sup, '%.17g');
else
    theta_sup_str = cell_data.theta_sup;
end

[area_bounds, perimeter_bounds] = compute_geometry_bounds(x_inf_str, x_sup_str, ...
    theta_inf_str, theta_sup_str);

fprintf('    Area bounds: [%.17f, %.17f]\n', area_bounds(1), area_bounds(2));
fprintf('    Perimeter bounds: [%.17f, %.17f]\n', perimeter_bounds(1), perimeter_bounds(2));

%% Step 3: Compute J lower bound
fprintf('  Computing J lower bound...\n');

[J_lower, J_diag] = compute_J_lower_bound(conjecture_type, lam1_lower(1), ...
    area_bounds, perimeter_bounds);

fprintf('    J lower bound: %.10e\n', J_lower);

%% Step 4: Determine verification status
verified = (J_lower > 0);

if verified
    fprintf('  Result: VERIFIED (J >= 0)\n');
else
    fprintf('  Result: NOT VERIFIED (J_lower = %.10e <= 0)\n', J_lower);
end

%% Prepare diagnostics
diagnostics = struct();
diagnostics.conjecture_type = conjecture_type;
diagnostics.cell_data = cell_data;
diagnostics.lam1_lower = lam1_lower(1);
diagnostics.area_bounds = area_bounds;
diagnostics.perimeter_bounds = perimeter_bounds;
diagnostics.J_lower = J_lower;
diagnostics.J_diagnostics = J_diag;
diagnostics.verified = verified;

end
