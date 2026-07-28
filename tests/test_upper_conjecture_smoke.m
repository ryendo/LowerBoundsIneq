function test_upper_conjecture_smoke()
%TEST_UPPER_CONJECTURE_SMOKE Minimal double and optional INTLAB smoke tests.
project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(project_root, 'src', 'upper_conjecture'));

thin = upper_verify_thin_sector('double', 0.06);
assert(strcmp(thin.rigor, 'exploratory_double'));
assert(thin.scaled_margin_lower > 0);

mesh = upper_make_reference_mesh(18);
matrices = upper_reference_matrices(mesh);
cell_def = struct('id', 1, 'x_lo', 0.70, 'x_hi', 0.702, ...
    'y_lo', 0.20, 'y_hi', 0.202);
r = upper_verify_trial_cell(cell_def, mesh, matrices, 'double');
assert(strcmp(r.rigor, 'exploratory_double'));
assert(strcmp(r.backend, 'P1_N18'));
assert(isfinite(r.margin_lower));
assert(r.margin_lower > 0);
fprintf('double smoke: thin_scaled_margin=%.6e cell_margin=%.6e\n', ...
    thin.scaled_margin_lower, r.margin_lower);

if strcmp(getenv('RUN_INTLAB_SMOKE'), '1')
    upper_prepare_runtime(project_root, 'interval');
    thin_i = upper_verify_thin_sector('interval', 0.06);
    assert(thin_i.verified);
    % Use a deliberately small cell with a comfortable positive margin.
    cell_i = struct('id', 2, 'x_lo', 0.70, 'x_hi', 0.7001, ...
        'y_lo', 0.20, 'y_hi', 0.2001);
    r_i = upper_verify_trial_cell(cell_i, mesh, matrices, 'interval');
    assert(r_i.verified);
    model4 = upper_make_high_order_model(4, 8, 'interval');
    upper_cell = struct('id', 3, 'x_lo', 0.5, 'x_hi', 0.502, ...
        'y_lo', 0.742, 'y_hi', 0.744);
    r4 = upper_verify_trial_cell(upper_cell, [], model4, 'interval');
    assert(r4.verified);
    assert(strcmp(r4.backend,'P4_N8'));
    fprintf(['interval smoke: thin_scaled_margin=%.6e ', ...
        'P1_margin=%.6e P4_upper_margin=%.6e\n'], ...
        thin_i.scaled_margin_lower, r_i.margin_lower, r4.margin_lower);
else
    fprintf(['interval smoke skipped; set RUN_INTLAB_SMOKE=1 and ', ...
        'INTLAB_ROOT=/path/to/Intlab_V12 to enable it.\n']);
end
end
