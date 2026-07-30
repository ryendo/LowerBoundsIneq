function test_muhat1_rayleigh_double()
%TEST_MUHAT1_RAYLEIGH_DOUBLE  Non-rigorous regression for the model quotient.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'src', 'degenerate'));

r = muhat1_rayleigh_cover( ...
    'mode', 'double', ...
    't', '0.38', ...
    'degree', 10, ...
    's_cells', 20);

assert(strcmp(r.rigor, 'exploratory_double'));
assert(r.max_midpoint_ritz > 13);
assert(r.max_midpoint_ritz < 15);
assert(~r.target_11_5_verified);
fprintf('double maximum midpoint Ritz value: %.12f\n', ...
    r.max_midpoint_ritz);
end
