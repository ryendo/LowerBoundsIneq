base_triangle = [0, 0, 1, 0, 0.5, sqrt(3)/2];

% Slightly perturbed triangle
triangle = [0, 0, 1, 0, 0.5, sqrt(3)/2];

% Perturbation direction (horizontal)
e_direction = [1, 0];

% Parameters
i = 1;
N_spectral = 2;
N_LG = 8;
N_rho = 32;
fem_ord_LG = 3;

fprintf('  Computing for perturbed triangle...\n');
[lami, dlami, ddlami_lb] = calc_ddlami_lower_bound(...
    i, base_triangle, triangle, e_direction, N_spectral, N_LG, N_rho, fem_ord_LG);

lami, dlami, ddlami_lb