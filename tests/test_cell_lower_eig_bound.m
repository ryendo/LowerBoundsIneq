%% TEST_CELL_LOWER_EIG_BOUND
% Minimal + robust test for cell_lower_eig_bound(cell_data)
% - Tests CR branch (isLG=0)
% - Optionally tests LG branch (isLG=1) if dependencies are available
%
% This test does NOT assume INTLAB is installed; it runs with INTERVAL_MODE=0.
% If your environment uses INTLAB wrappers (I_intval/I_mid/I_inf), they must be on path.

clear; clc;

fprintf('================================================================\n');
fprintf('TEST: cell_lower_eig_bound (CR/LG lower eigenvalue bounds)\n');
fprintf('================================================================\n\n');

%% Paths (adjust to your repo layout)
% These should include:
% - cell_lower_eig_bound.m
% - I_intval, I_mid, I_inf (and INTLAB helpers if used)
% - make_mesh_by_gmsh, create_matrix_crouzeix_raviart, find_tri2edge, find_mesh_hmax
% - (optional LG) laplace_eig_lagrange_detailed, RT_Hdiv_problem_dirichlet, I_eig, I_hull, veigs
%
% Example:
% addpath('../src');
% addpath('../src/algorithms');
% addpath('../src/mesh');
% addpath('../src/fem');
% addpath('../src/interval');
% addpath('../src/lib/eigenvalue_bound');
% addpath('../src/lib/mesh');

%% Force non-interval mode (standard floating point)
global INTERVAL_MODE;
INTERVAL_MODE = 0;
fprintf('[Test] INTERVAL_MODE = %d\n\n', INTERVAL_MODE);

%% Helper: build a "cell_data" struct
% We choose a geometry cell around the equilateral triangle:
% x ~ 0.5, theta ~ pi/3  (since y = x * tan(theta), so y= sqrt(3)/2 at x=0.5, theta=pi/3)
%
% IMPORTANT: x_inf/x_sup/theta_inf/theta_sup are strings in your pipeline.
cell_data = struct();
cell_data.x_inf       = '0.5';
cell_data.x_sup       = '0.5';
cell_data.theta_inf   = '1.0471975511965976';  % pi/3
cell_data.theta_sup   = '1.0471975511965976';  % pi/3

% Mesh/FEM params (CR)
cell_data.mesh_size_lower_cr = 0.25;  % coarse for speed

% --- CR branch ---
cell_data.isLG = 0;

% Request only lambda_1 by default; also test neig=2 as a second case
cell_data.neig = 1;

%% -----------------------------
% Test 1: CR branch, neig=1
% -----------------------------
fprintf('----------------------------------------------------------------\n');
fprintf('[Test 1] CR branch (isLG=0), neig=1\n');
fprintf('----------------------------------------------------------------\n');

test1_passed = false;

try
    eig_bounds = cell_lower_eig_bound(cell_data);

    % Basic shape checks
    assert(isnumeric(eig_bounds) && isvector(eig_bounds), 'eig_bounds must be a numeric vector.');
    assert(numel(eig_bounds) == 1, 'Expected 1 eigenvalue bound (neig=1).');
    assert(isfinite(eig_bounds(1)), 'eig_bounds(1) must be finite.');
    assert(eig_bounds(1) > 0, 'Lower bound must be positive.');

    % Weak sanity check vs known equilateral lambda_1 = 16*pi^2/3 (~52.64)
    lambda1_exact = 16*pi^2/3;
    fprintf('[Test 1] lambda1_exact (equilateral) = %.6f\n', lambda1_exact);
    fprintf('[Test 1] computed lower bound        = %.6f\n', eig_bounds(1));

    % A validated lower bound should be <= exact (for same domain).
    % In practice, if the geometry is exactly equilateral and the method is correct,
    % the CR validated lower bound should not exceed the true eigenvalue.
    assert(eig_bounds(1) <= lambda1_exact * 1.05, ...
        'Lower bound looks too large (unexpected for a lower bound).');

    test1_passed = true;
    fprintf('[Test 1] ✓ PASSED\n\n');

catch ME
    fprintf('[Test 1] ✗ FAILED: %s\n', ME.message);
    for k = 1:numel(ME.stack)
        fprintf('    %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
    end
    fprintf('\n');
end

%% -----------------------------
% Test 2: CR branch, neig=2
% -----------------------------
fprintf('----------------------------------------------------------------\n');
fprintf('[Test 2] CR branch (isLG=0), neig=2\n');
fprintf('----------------------------------------------------------------\n');

test2_passed = false;

try
    cell_data2 = cell_data;
    cell_data2.neig = 2;

    eig_bounds2 = cell_lower_eig_bound(cell_data2);

    assert(isnumeric(eig_bounds2) && isvector(eig_bounds2), 'eig_bounds must be a numeric vector.');
    assert(numel(eig_bounds2) == 2, 'Expected 2 eigenvalue bounds (neig=2).');
    assert(all(isfinite(eig_bounds2)), 'All bounds must be finite.');
    assert(all(eig_bounds2 > 0), 'All bounds must be positive.');
    assert(eig_bounds2(2) >= eig_bounds2(1), 'Eigenvalue bounds should be nondecreasing.');

    fprintf('[Test 2] lower bounds: [%.6f, %.6f]\n', eig_bounds2(1), eig_bounds2(2));
    test2_passed = true;
    fprintf('[Test 2] ✓ PASSED\n\n');

catch ME
    fprintf('[Test 2] ✗ FAILED: %s\n', ME.message);
    for k = 1:numel(ME.stack)
        fprintf('    %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
    end
    fprintf('\n');
end

%% -----------------------------
% Test 3 (optional): LG branch smoke test
% -----------------------------
% LG needs many dependencies and is heavier. Keep it optional.
RUN_LG = true;

test3_passed = false; % remains true if skipped

if RUN_LG
    fprintf('----------------------------------------------------------------\n');
    fprintf('[Test 3] LG branch (isLG=1) SMOKE TEST\n');
    fprintf('----------------------------------------------------------------\n');

    try
        cell_dataLG = cell_data;
        cell_dataLG.isLG = 1;
        cell_dataLG.neig = 1;

        % LG params required by your function
        cell_dataLG.mesh_size_lower_LG  = 0.4;   % coarse for speed
        cell_dataLG.fem_order_lower_LG  = '3';   % string as expected by str2num

        eig_boundsLG = cell_lower_eig_bound(cell_dataLG);

        assert(isnumeric(eig_boundsLG) && isvector(eig_boundsLG), 'eig_bounds must be a numeric vector.');
        assert(numel(eig_boundsLG) == 1, 'Expected 1 eigenvalue bound (neig=1).');
        assert(isfinite(eig_boundsLG(1)) && eig_boundsLG(1) > 0, 'LG bound must be positive and finite.');

        fprintf('[Test 3] LG lower bound = %.6f\n', eig_boundsLG(1));
        test3_passed = true;
        fprintf('[Test 3] ✓ PASSED\n\n');

    catch ME
        fprintf('[Test 3] ✗ FAILED: %s\n', ME.message);
        for k = 1:numel(ME.stack)
            fprintf('    %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
        end
        fprintf('\n');
        test3_passed = false;
    end
else
    fprintf('----------------------------------------------------------------\n');
    fprintf('[Test 3] LG branch skipped (set RUN_LG=true to enable)\n');
    fprintf('----------------------------------------------------------------\n\n');
end

%% Summary
fprintf('================================================================\n');
fprintf('TEST SUMMARY\n');
fprintf('================================================================\n');
fprintf('Test 1 (CR, neig=1): %s\n', ternary(test1_passed, 'PASSED', 'FAILED'));
fprintf('Test 2 (CR, neig=2): %s\n', ternary(test2_passed, 'PASSED', 'FAILED'));
fprintf('Test 3 (LG smoke):   %s\n', ternary(test3_passed, 'PASSED', 'FAILED'));
fprintf('----------------------------------------------------------------\n');

all_ok = test1_passed && test2_passed && test3_passed;
if all_ok
    fprintf('✓✓✓ ALL TESTS PASSED ✓✓✓\n');
else
    fprintf('✗✗✗ SOME TESTS FAILED ✗✗✗\n');
end
fprintf('================================================================\n\n');

%% Local helper
function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end
