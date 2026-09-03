%RUN_VER10_RELEASE_TESTS  Consolidated release regression suite.

project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(project_root,'tests'));
addpath(fullfile(project_root,'scripts_run'));
if ~strcmp(getenv('RUN_INTLAB_SMOKE'),'1')
    error('run_ver10_release_tests:IntlabSmokeRequired', ...
        ['Set RUN_INTLAB_SMOKE=1.  The release suite must not report ', ...
         'success after skipping the interval upper-bound tests.']);
end
omega_up_all_prepare_worker(project_root,'interval');

test_verified_ritz_enclosures(true);
test_cell_lower_eig_bound;
test_residual_hessian_estimator_double(true);
test_residual_hessian_general_index;
test_bernstein_strong_residual_estimator(true);
test_omega_up_spectral_atlas;
test_omega_up_all_functionals_unified;
test_omega_up_checkpoint_validation;
test_omega_mid_geometric_coverage;
test_omega_mid_artifact_binding;
test_omega_mid_floor_geometry;
test_omega_mid_unified_driver;
test_apply_exact_boundary_point_setting;
test_upper_conjecture_smoke;
test_muhat1_rayleigh_certificate;
test_upper_compact_cover;
test_upper_global_manifest_combiner;
% This test deliberately replaces and deletes the current parallel pool,
% so keep it last in the consolidated suite.
test_worker_runtime_binding;

fprintf('run_ver10_release_tests: ALL PASS\n');
