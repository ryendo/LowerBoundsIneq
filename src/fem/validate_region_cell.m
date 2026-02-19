function cell_result = validate_region_cell(cell_data, varargin)
% VALIDATE_REGION_CELL: Validate a single cell for eigenvalue/functional bounds
%
% Paper Reference: Section 4, Step 2 (Omega_mid verification)
%   This function supports two modes:
%   1. J evaluation mode (default): Verify J1 or J2 >= 0
%   2. Gap check mode (legacy): Verify lambda_2 < lambda_3
%
% Inputs:
%   cell_data: structure with geometry and FEM parameters
%   varargin: optional parameters
%     'mode': 'J' (functional) or 'gap' (eigenvalue gap), default: 'J'
%     'conjecture_type': 'J1' or 'J2' (only for mode='J'), default: 'J1'
%
% Outputs:
%   cell_result: structure with validation results
%     For mode='J':
%       - J_lower: lower bound on J
%       - verified: true if J_lower > 0
%     For mode='gap':
%       - lam2_sup, lam3_inf: eigenvalue bounds
%       - verified: true if lam3_inf > lam2_sup
%
% Author: Based on paper by R. Endo, X. Liu, P. Mariano
% Date: 2025-01-14 (modified for J1/J2 evaluation)

% Parse optional parameters
p = inputParser;
addParameter(p, 'mode', 'J', @ischar);
addParameter(p, 'conjecture_type', 'J1', @ischar);
parse(p, varargin{:});

mode = p.Results.mode;
conjecture_type = p.Results.conjecture_type;

if strcmpi(mode, 'J')
    %% ====================================================================
    %  J Evaluation Mode (New - for J1/J2 functional)
    %  ====================================================================

    % Compute J lower bound using verify_J_positive
    % Note: This requires Paper_Algorithms to be on the path
    addpath('../Paper_Algorithms');

    fprintf('Validating cell (J mode, %s)...\n', conjecture_type);

    [verified, J_lower, diagnostics] = verify_J_positive(conjecture_type, cell_data);

    % Construct result structure
    cell_result = struct();
    cell_result.i = cell_data.i;
    cell_result.x_inf = cell_data.x_inf;
    cell_result.x_sup = cell_data.x_sup;
    cell_result.theta_inf = cell_data.theta_inf;
    cell_result.theta_sup = cell_data.theta_sup;
    cell_result.mode = 'J';
    cell_result.conjecture_type = conjecture_type;
    cell_result.J_lower = J_lower;
    cell_result.verified = verified;
    cell_result.diagnostics = diagnostics;

    if verified
        fprintf('Result: OK (J >= 0, J_lower = %.6e)\n', J_lower);
    else
        fprintf('Result: NG (J_lower = %.6e <= 0)\n', J_lower);
    end

elseif strcmpi(mode, 'gap')
    %% ====================================================================
    %  Gap Check Mode (Legacy - for lambda_2 < lambda_3)
    %  ====================================================================

    fprintf('Validating cell (gap mode)...\n');

    % Set neig = 3 for gap check
    cell_data_gap = cell_data;
    cell_data_gap.neig = 3;

    % Compute Upper Bound
    tic;
    disp('  Upper bound computation...');
    cell_ub = cell_upper_eig_bound(cell_data_gap);
    toc;

    % Compute Lower Bound
    tic;
    disp('  Lower bound computation...');
    cell_lb = cell_lower_eig_bound(cell_data_gap);
    toc;

    % Validation Logic (Check if Lower(3) > Upper(2))
    disp('  Region cell validation:');

    if cell_lb(3) > cell_ub(2)
        disp('    OK');
        verified = true;
    else
        disp('    NG');
        fprintf('    Gap Closed! LB(3)=%.17g, UB(2)=%.17g\n', cell_lb(3), cell_ub(2));
        verified = false;
    end

    % Construct the result structure
    cell_result = struct();
    cell_result.i = cell_data.i;
    cell_result.x_inf = cell_data.x_inf;
    cell_result.x_sup = cell_data.x_sup;
    cell_result.theta_inf = cell_data.theta_inf;
    cell_result.theta_sup = cell_data.theta_sup;
    cell_result.mode = 'gap';

    % Format numerical results to high-precision strings
    cell_result.lam2_sup = compose("%.17g", cell_ub(2));
    cell_result.lam3_inf = compose("%.17g", cell_lb(3));
    cell_result.verified = verified;

else
    error('Unknown mode: %s. Use ''J'' or ''gap''.', mode);
end

end
