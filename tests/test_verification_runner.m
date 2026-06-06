%% TEST_VERIFICATION_RUNNER
% Test the VerificationRunner class

clear; clc;

fprintf('================================================================\n');
fprintf('TEST: VerificationRunner Class\n');
fprintf('================================================================\n\n');

%% Test 1: Initialization
fprintf('[Test 1] Testing initialization...\n');
try
    runner = VerificationRunner('verbose', true, 'mesh_size', 0.1, 'grid_resolution', 5);
    fprintf('[Test 1] ✓ Initialization successful\n\n');
    test1_passed = true;
catch ME
    fprintf('[Test 1] ✗ Initialization failed: %s\n\n', ME.message);
    test1_passed = false;
end

if ~test1_passed
    fprintf('Cannot continue without successful initialization\n');
    return;
end

%% Test 2: Compute eigenvalue bounds
fprintf('[Test 2] Testing eigenvalue computation...\n');
try
    triangle_eq = [0, 0, 1, 0, 0.5, sqrt(3)/2];
    [lambda_bounds, diag] = runner.computeEigenvalueBounds(triangle_eq, 5);

    fprintf('[Test 2] ✓ Eigenvalue computation successful\n');
    fprintf('[Test 2] Computed %d eigenvalues\n', length(lambda_bounds));
    fprintf('[Test 2] λ_1 ≈ %.6f (exact: %.6f)\n', lambda_bounds(1), 16*pi^2/3);
    fprintf('[Test 2] Relative error: %.2e\n\n', ...
        abs(lambda_bounds(1) - 16*pi^2/3) / (16*pi^2/3));

    test2_passed = true;
catch ME
    fprintf('[Test 2] ✗ Eigenvalue computation failed: %s\n\n', ME.message);
    test2_passed = false;
end

%% Test 3: Verify Omega_up
fprintf('[Test 3] Testing Omega_up verification...\n');
try
    results_up = runner.verifyOmegaUp('J1');
    fprintf('[Test 3] ✓ Omega_up verification completed\n');
    fprintf('[Test 3] Result: %s\n\n', mat2str(results_up.verified));
    test3_passed = true;
catch ME
    fprintf('[Test 3] ✗ Omega_up verification failed: %s\n\n', ME.message);
    test3_passed = false;
end

%% Test 4: Verify Omega_mid (quick)
fprintf('[Test 4] Testing Omega_mid verification...\n');
try
    results_mid = runner.verifyOmegaMid('J1');
    fprintf('[Test 4] ✓ Omega_mid verification completed\n');
    fprintf('[Test 4] Grid points: %d/%d verified\n', ...
        results_mid.n_verified, results_mid.n_total);
    fprintf('[Test 4] Result: %s\n\n', mat2str(results_mid.verified));
    test4_passed = true;
catch ME
    fprintf('[Test 4] ✗ Omega_mid verification failed: %s\n\n', ME.message);
    test4_passed = false;
end

%% Test 5: Verify Omega_down
fprintf('[Test 5] Testing Omega_down verification...\n');
try
    results_down = runner.verifyOmegaDown('J1');
    fprintf('[Test 5] ✓ Omega_down verification completed\n');
    fprintf('[Test 5] Result: %s\n\n', mat2str(results_down.verified));
    test5_passed = true;
catch ME
    fprintf('[Test 5] ✗ Omega_down verification failed: %s\n\n', ME.message);
    test5_passed = false;
end

%% Test 6: Complete verification
fprintf('[Test 6] Testing complete verification workflow...\n');
try
    results_complete = runner.runCompleteVerification('J1');
    fprintf('[Test 6] ✓ Complete verification finished\n');
    fprintf('[Test 6] Overall result: %s\n', mat2str(results_complete.overall_verified));
    fprintf('[Test 6] Total time: %.2f seconds\n\n', results_complete.total_time);
    test6_passed = true;
catch ME
    fprintf('[Test 6] ✗ Complete verification failed: %s\n\n', ME.message);
    test6_passed = false;
end

%% Test 7: Report generation
fprintf('[Test 7] Testing report generation...\n');
try
    runner.generateReport();
    fprintf('[Test 7] ✓ Report generated successfully\n\n');
    test7_passed = true;
catch ME
    fprintf('[Test 7] ✗ Report generation failed: %s\n\n', ME.message);
    test7_passed = false;
end

%% Summary
fprintf('================================================================\n');
fprintf('TEST SUMMARY\n');
fprintf('================================================================\n');
fprintf('Test 1 (Initialization):        %s\n', ternary(test1_passed, '✓ PASSED', '✗ FAILED'));
fprintf('Test 2 (Eigenvalue bounds):     %s\n', ternary(test2_passed, '✓ PASSED', '✗ FAILED'));
fprintf('Test 3 (Omega_up):              %s\n', ternary(test3_passed, '✓ PASSED', '✗ FAILED'));
fprintf('Test 4 (Omega_mid):             %s\n', ternary(test4_passed, '✓ PASSED', '✗ FAILED'));
fprintf('Test 5 (Omega_down):            %s\n', ternary(test5_passed, '✓ PASSED', '✗ FAILED'));
fprintf('Test 6 (Complete verification): %s\n', ternary(test6_passed, '✓ PASSED', '✗ FAILED'));
fprintf('Test 7 (Report generation):     %s\n', ternary(test7_passed, '✓ PASSED', '✗ FAILED'));

all_passed = test1_passed && test2_passed && test3_passed && test4_passed && ...
             test5_passed && test6_passed && test7_passed;

fprintf('================================================================\n');
if all_passed
    fprintf('✓✓✓ ALL TESTS PASSED ✓✓✓\n');
else
    fprintf('✗✗✗ SOME TESTS FAILED ✗✗✗\n');
end
fprintf('================================================================\n\n');

%% Helper function
function result = ternary(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end
