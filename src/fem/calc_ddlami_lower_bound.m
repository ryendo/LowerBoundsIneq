function [lami, dlami, ddlam_i_lower_bound, diagnostics] = calc_ddlami_lower_bound(i, base_triangle, triangle, e_direction, N_spectral, N_LG, N_rho, fem_ord_LG, M_neumann)
%CALC_DDLAMI_LOWER_BOUND  Certified lower bound for ddot(lambda_i).
%
% Paper reference:
%   Theorem "Certified FEM Dirichlet--Neumann estimator"
%   (label thm:main-theorem-est).
%
% This routine evaluates the lower side of the paper estimate
%
%   underline D_{i,N}^{p,h}
%   - 2 * lambda_{N+1}^p/(lambda_{N+1}^p-lambda_i^p)
%       * (overline U_{i,M}^{p,h} - underline S_{i,N}^{p,h})
%   <= ddot(lambda_i^p).
%
% The finite Dirichlet quantity D_{i,N} is not used alone as a lower bound.
% Its missing negative tail is bounded by the Neumann complement U_{i,M};
% see equations defining D, S, U and Theorem main-theorem-est in the paper.

if nargin < 9 || isempty(M_neumann)
    M_neumann = max(N_spectral, 2);
end

if i > N_spectral
    error('The paper bound requires i <= N_spectral.');
end
if i ~= 1
    error('The eigenspace-based paper estimator is implemented for lambda_1.');
end
if M_neumann < 1
    error('M_neumann must be at least 1.');
end

% We need lambda_{N+1} for the tail coefficient and one extra eigenvalue for
% the Liu--Vejchodsky eigenspace estimator when the last cluster is formed.
num_dirichlet_needed = max(N_spectral + 3, N_spectral + 1);

% -------------------------------------------------------------------------
% Dirichlet FEM data: eigenvalue enclosures, conforming eigenvectors, and
% eta_phi in Theorem main-theorem-est.
% -------------------------------------------------------------------------
[lams_raw, uh_list, K_D, M_D, A_xx, A_xy, A_yy, uh_list_full, meshCG] = ...
    calc_eigen_bounds_any_order_1k_wh(num_dirichlet_needed, base_triangle, N_LG, N_rho, fem_ord_LG);

lams_h_raw = I_intval(zeros(num_dirichlet_needed, 1));
for k = 1:num_dirichlet_needed
    lams_h_raw(k) = (uh_list(:,k)' * K_D * uh_list(:,k)) / (uh_list(:,k)' * M_D * uh_list(:,k));
end

x_base   = I_intval(base_triangle(5));
y_base   = I_intval(base_triangle(6));
x_target = I_intval(triangle(5));
y_target = I_intval(triangle(6));

% Eigenvalue perturbation over the target cell, used exactly as in the
% paper's uniform cell algorithms for Omega_up.
S_inv     = [1, (x_base - x_target) / y_target; 0, y_base / y_target];
[eig_factor_lower, eig_factor_upper] = affine_metric_eigenvalue_factors(S_inv);

lams_interval = scale_positive_eigen_bounds(lams_raw, eig_factor_lower, eig_factor_upper);
lams_h        = scale_positive_eigen_bounds(lams_h_raw, eig_factor_lower, eig_factor_upper);

dirichlet_clusters = auto_cluster_eigenvalues(lams_interval, 0.01);
[eta_phi_subspace, eta_phi_cluster_sq, eta_phi_cluster_b_sq, eta_phi_info] = ...
    calc_grad_error_bounds(lams_interval, lams_h, uh_list, K_D, M_D, dirichlet_clusters);
dirichlet_cluster_ids = complete_cluster_ids(dirichlet_clusters, N_spectral, 'Dirichlet');
if ~isequal(dirichlet_clusters{dirichlet_cluster_ids(1)}, 1)
    error('The eigenspace estimator requires lambda_1 to be an isolated first cluster.');
end
eta_phi_1 = individual_est_a_from_singleton( ...
    lams_interval(1), lams_h(1), eta_phi_cluster_b_sq(dirichlet_cluster_ids(1)));

% -------------------------------------------------------------------------
% Neumann FEM data: positive modes mu_1,...,mu_M plus one extra mode for the
% eigenspace estimator eta_psi.  These are used only in U_{i,M}.
% -------------------------------------------------------------------------
num_neumann_needed = M_neumann + 1;
[mu_raw, mu_h_raw, psi_list, K_N, M_N, N_xx, N_xy, N_yy, N_ux_vy, N_uy_vx] = ...
    calc_neumann_eigen_bounds_any_order_1k_wh(num_neumann_needed, base_triangle, N_LG, N_rho, fem_ord_LG, meshCG);

if size(uh_list_full, 1) ~= size(psi_list, 1)
    error('Dirichlet full-dof vectors and Neumann vectors live on different meshes.');
end

% The perturbation estimate applies to both Dirichlet and Neumann eigenvalues
% under the same affine cell map; see the paper's eigenvalue perturbation
% lemma for either boundary condition.  The Neumann eigenspace errors used in
% epsilon_C are therefore computed from the scaled cell enclosures.
mu_interval = scale_positive_eigen_bounds(mu_raw, eig_factor_lower, eig_factor_upper);
mu_h = scale_positive_eigen_bounds(mu_h_raw, eig_factor_lower, eig_factor_upper);
mu_clusters = auto_cluster_eigenvalues(mu_interval, 0.01);
[eta_psi, eta_psi_cluster_sq, eta_psi_cluster_b_sq, eta_psi_info] = ...
    calc_grad_error_bounds(mu_interval, mu_h, psi_list, K_N, M_N, mu_clusters);

% -------------------------------------------------------------------------
% Shape derivative matrices from the paper's first/second derivative
% formulas for vertex perturbations.
% -------------------------------------------------------------------------
a = e_direction(1);
b = e_direction(2);

dotP  = [[0,                    -a/y_target]; ...
         [-a/y_target,          -2*b/y_target]];
ddotP = [[2*a^2/y_target^2,      4*a*b/y_target^2]; ...
         [4*a*b/y_target^2,      6*b^2/y_target^2]];

norm_dotP  = norm(dotP, 2);
norm_ddotP = norm(ddotP, 2);

shape_D = shape_matrix(dotP, A_xx, A_xy, A_yy);
shape_D2 = shape_matrix(ddotP, A_xx, A_xy, A_yy);

lami  = lams_interval(i);
lamih = lams_h(i);

% -------------------------------------------------------------------------
% First derivative enclosure, used by Omega_up y-direction via J_yy-simple.
% This is Lemma "A priori error estimate for the discrete first-order shape
% derivative".
% -------------------------------------------------------------------------
dlamih = uh_list(:,i)' * shape_D * uh_list(:,i);
e_dlami = norm_dotP * (sqrt(lami) + sqrt(lamih)) * eta_phi_1;
dlami = I_hull(dlamih - e_dlami, dlamih + e_dlami);

% -------------------------------------------------------------------------
% A_i and epsilon_A: equations explicit-eps-A and explicit-A-bounds.
% -------------------------------------------------------------------------
Ahat_i = uh_list(:,i)' * shape_D2 * uh_list(:,i);
eps_A = norm_ddotP * eta_phi_1 * (sqrt(lami) + sqrt(lamih));
A_lower = Ahat_i - eps_A;
A_upper = Ahat_i + eps_A;

% -------------------------------------------------------------------------
% Eigenspace-based Dirichlet quantities in the revised paper:
%   equations D-lower-eigenspace, D-upper-eigenspace, and
%   S-lower-eigenspace.
%
% Each complete consecutive Dirichlet range E_k is evaluated through the
% basis-invariant projection norm
%   || P_{grad E_k^h}(dotP grad phihat_1) ||.
% This avoids labeling eigenfunctions inside a multiple/clustered range.
% -------------------------------------------------------------------------
eta_G = norm_dotP * eta_phi_1;
Ghat = gradient_image_norm_sq(uh_list(:,i), dotP, A_xx, A_xy, A_yy);
eps_G = norm_dotP^2 * eta_phi_1 * (sqrt(lami) + sqrt(lamih));
G_upper = Ghat + eps_G;

D_lower = A_lower;
D_upper = A_upper;
S_lower = I_intval(0);

num_dirichlet_blocks = numel(dirichlet_cluster_ids);
dir_proj_lower = I_intval(zeros(num_dirichlet_blocks, 1));
dir_proj_upper = I_intval(zeros(num_dirichlet_blocks, 1));
dir_proj_radius = I_intval(zeros(num_dirichlet_blocks, 1));
dir_proj_energy_lower = I_intval(zeros(num_dirichlet_blocks, 1));
dir_proj_energy_upper = I_intval(zeros(num_dirichlet_blocks, 1));
dir_proj_coeff = cell(num_dirichlet_blocks, 1);

for q = 1:num_dirichlet_blocks
    cluster_id = dirichlet_cluster_ids(q);
    idx_cluster = dirichlet_clusters{cluster_id};
    idx_cluster = idx_cluster(:).';

    [proj_lower, proj_upper, coeff] = projection_norm_bounds( ...
        uh_list(:,i), uh_list(:,idx_cluster), shape_D, K_D);

    radius = sqrt(eta_phi_cluster_sq(cluster_id)) * sqrt(G_upper) + eta_G;
    energy_lower = squared_positive_part_lower(proj_lower, radius);
    energy_upper = squared_sum_upper(proj_upper, radius);

    dir_proj_lower(q) = proj_lower;
    dir_proj_upper(q) = proj_upper;
    dir_proj_radius(q) = radius;
    dir_proj_energy_lower(q) = energy_lower;
    dir_proj_energy_upper(q) = energy_upper;
    dir_proj_coeff{q} = coeff;

    S_lower = S_lower + energy_lower;

    if isequal(idx_cluster, 1)
        continue;
    end

    lam_n_lower = I_inf(lams_interval(idx_cluster(1)));
    lam_N_upper = I_sup(lams_interval(idx_cluster(end)));
    lam_1_lower = I_inf(lams_interval(1));
    lam_1_upper = I_sup(lams_interval(1));

    if ~(lam_n_lower > lam_1_upper)
        error('Cannot certify lambda_%d > lambda_1 for a Dirichlet cluster factor.', idx_cluster(1));
    end
    if ~(lam_N_upper > lam_1_lower)
        error('Cannot certify a positive Dirichlet cluster denominator.');
    end

    factor_upper = I_intval(lam_n_lower / (lam_n_lower - lam_1_upper));
    factor_lower = I_intval(lam_N_upper / (lam_N_upper - lam_1_lower));
    D_lower = D_lower - 2 * factor_upper * energy_upper;
    D_upper = D_upper - 2 * factor_lower * energy_lower;
end

% -------------------------------------------------------------------------
% Eigenspace-based Neumann complement: equation U-upper-eigenspace.
%
%   || P_{J grad E_K^N} (dotP grad phi_i) ||^2
%     = sum_{m in K} |C_im|^2/mu_m.
%
% The revised paper evaluates this projection norm directly for each complete
% Neumann range F_l; no individual Neumann eigenfunction labels are used.
% -------------------------------------------------------------------------
U_upper = G_upper;
C_hat = I_intval(zeros(M_neumann, 1));
neumann_coupling = neumann_complement_matrix(dotP, N_xx, N_yy, N_ux_vy, N_uy_vx);
for m = 1:M_neumann
    C_hat(m) = uh_list_full(:,i)' * neumann_coupling * psi_list(:,m);
end

neumann_cluster_ids = complete_cluster_ids(mu_clusters, M_neumann, 'Neumann');
num_neumann_blocks = numel(neumann_cluster_ids);
neu_proj_lower = I_intval(zeros(num_neumann_blocks, 1));
neu_proj_upper = I_intval(zeros(num_neumann_blocks, 1));
neu_proj_radius = I_intval(zeros(num_neumann_blocks, 1));
neu_proj_energy_lower = I_intval(zeros(num_neumann_blocks, 1));
neu_proj_coeff = cell(num_neumann_blocks, 1);

for q = 1:num_neumann_blocks
    cluster_id = neumann_cluster_ids(q);
    idx_cluster = mu_clusters{cluster_id};
    idx_cluster = idx_cluster(:).';

    if any(I_inf(mu_interval(idx_cluster)) <= 0)
        error('A Neumann cluster used in U contains a non-positive eigenvalue lower bound.');
    end

    [proj_lower, proj_upper, coeff] = projection_norm_bounds( ...
        uh_list_full(:,i), psi_list(:,idx_cluster), neumann_coupling, K_N);

    radius = sqrt(eta_psi_cluster_sq(cluster_id)) * sqrt(G_upper) + eta_G;
    energy_lower = squared_positive_part_lower(proj_lower, radius);

    neu_proj_lower(q) = proj_lower;
    neu_proj_upper(q) = proj_upper;
    neu_proj_radius(q) = radius;
    neu_proj_energy_lower(q) = energy_lower;
    neu_proj_coeff{q} = coeff;

    U_upper = U_upper - energy_lower;
end

% -------------------------------------------------------------------------
% Tail coefficient lambda_{N+1}/(lambda_{N+1}-lambda_i).  The paper's
% Omega_up algorithm uses lower(lambda_{N+1}) and upper(lambda_i) to obtain
% a certified upper bound for this decreasing coefficient.
% -------------------------------------------------------------------------
lamNp1_lower = I_inf(lams_interval(N_spectral + 1));
lami_upper = I_sup(lams_interval(i));
if ~(lamNp1_lower > lami_upper)
    error('Cannot certify lambda_{N+1} > lambda_i for the Neumann tail factor.');
end
Q_tail_upper = I_intval(lamNp1_lower / (lamNp1_lower - lami_upper));

tail_gap_upper = I_sup(U_upper) - I_inf(S_lower);
if tail_gap_upper < 0
    error('Computed U upper bound is below S lower bound; cannot certify the DN tail.');
end

ddlam_i_lower_bound = I_intval(I_inf(D_lower - 2 * Q_tail_upper * I_intval(tail_gap_upper)));

diagnostics = struct();
diagnostics.N_spectral = N_spectral;
diagnostics.M_neumann = M_neumann;
diagnostics.lami = lami;
diagnostics.dlami = dlami;
diagnostics.eig_factor_lower = eig_factor_lower;
diagnostics.eig_factor_upper = eig_factor_upper;
diagnostics.D_lower = D_lower;
diagnostics.D_upper = D_upper;
diagnostics.ddlam_lower = ddlam_i_lower_bound;
diagnostics.ddlam_upper = I_intval(I_sup(D_upper));
diagnostics.ddlam_width = I_sup(D_upper) - I_inf(ddlam_i_lower_bound);
diagnostics.S_lower = S_lower;
diagnostics.U_upper = U_upper;
diagnostics.Ghat = Ghat;
diagnostics.G_upper = G_upper;
diagnostics.Q_tail_upper = Q_tail_upper;
diagnostics.tail_gap_upper = tail_gap_upper;
diagnostics.eta_phi_i = eta_phi_1;
diagnostics.eta_phi_1 = eta_phi_1;
diagnostics.eta_phi_subspace = eta_phi_subspace(1:max(N_spectral, 1));
diagnostics.eta_phi_cluster = sqrt(eta_phi_cluster_sq);
diagnostics.eta_phi_info = eta_phi_info;
diagnostics.dirichlet_clusters = dirichlet_clusters;
diagnostics.dirichlet_cluster_ids = dirichlet_cluster_ids;
diagnostics.eta_G = eta_G;
diagnostics.dir_proj_lower = dir_proj_lower;
diagnostics.dir_proj_upper = dir_proj_upper;
diagnostics.dir_proj_radius = dir_proj_radius;
diagnostics.dir_proj_energy_lower = dir_proj_energy_lower;
diagnostics.dir_proj_energy_upper = dir_proj_energy_upper;
diagnostics.dir_proj_coeff = dir_proj_coeff;
diagnostics.eta_psi = sqrt(eta_psi_cluster_sq);
diagnostics.eta_psi_cluster = sqrt(eta_psi_cluster_sq);
diagnostics.eta_psi_cluster_sq = eta_psi_cluster_sq;
diagnostics.eta_psi_cluster_b = sqrt(eta_psi_cluster_b_sq);
diagnostics.eta_psi_cluster_b_sq = eta_psi_cluster_b_sq;
diagnostics.eta_psi_info = eta_psi_info;
diagnostics.mu_clusters = mu_clusters;
diagnostics.neumann_cluster_ids = neumann_cluster_ids;
diagnostics.mu_interval = mu_interval(1:M_neumann);
diagnostics.mu_h = mu_h(1:M_neumann);
diagnostics.C_hat = C_hat;
diagnostics.neu_proj_lower = neu_proj_lower;
diagnostics.neu_proj_upper = neu_proj_upper;
diagnostics.neu_proj_radius = neu_proj_radius;
diagnostics.neu_proj_energy_lower = neu_proj_energy_lower;
diagnostics.neu_proj_coeff = neu_proj_coeff;

end


function A = shape_matrix(P, A_xx, A_xy, A_yy)
% Matrix for (P grad u, grad v).
A = P(1,1)*A_xx + (P(1,2)+P(2,1))*A_xy + P(2,2)*A_yy;
end


function val = gradient_image_norm_sq(u, P, A_xx, A_xy, A_yy)
% Discrete Ghat = ||P grad u||^2, equation G-Ghat-fem.
c_xx = P(1,1)^2 + P(2,1)^2;
c_xy = P(1,1)*P(1,2) + P(2,1)*P(2,2);
c_yy = P(1,2)^2 + P(2,2)^2;
val = u' * (c_xx*A_xx + 2*c_xy*A_xy + c_yy*A_yy) * u;
end


function [factor_lower, factor_upper] = affine_metric_eigenvalue_factors(S_inv)
% Eigenvalue perturbation factors for the affine cell map.
%
% Lemma eig-perturbation-first gives
%   lambda_min(S^{-1}S^{-T}) * nu(base)
%     <= nu(target)
%     <= lambda_max(S^{-1}S^{-T}) * nu(base).
% The lower and upper factors must be kept separate; using only the spectral
% norm makes the cell enclosure look artificially narrow and can corrupt the
% Liu--Vejchodsky eigenspace error estimates.
%
% In our cell algorithms S_inv has the triangular form [1,a;0,b], with
% a containing 0 and b > 0.  Direct interval eig formulas suffer severe
% dependency overestimation, so evaluate the scalar 2x2 eigenvalue formula at
% the endpoint candidates for |a| and b.
aI = S_inv(1,2);
bI = S_inv(2,2);

alpha_max = max(abs(I_inf(aI)), abs(I_sup(aI)));
beta_min = min(abs(I_inf(bI)), abs(I_sup(bI)));
beta_max = max(abs(I_inf(bI)), abs(I_sup(bI)));

lambda_min_candidates = [];
lambda_max_candidates = [];
for alpha = [0, alpha_max]
    for beta = [beta_min, beta_max]
        [emin, emax] = triangular_metric_eigs(alpha, beta);
        lambda_min_candidates(end+1) = emin; %#ok<AGROW>
        lambda_max_candidates(end+1) = emax; %#ok<AGROW>
    end
end

factor_lower = I_intval(max(min(lambda_min_candidates), 0));
factor_upper = I_intval(max(lambda_max_candidates));
end


function [emin, emax] = triangular_metric_eigs(alpha, beta)
% Eigenvalues of [1,alpha;0,beta] * [1,alpha;0,beta]^T.
a11 = 1 + alpha^2;
a12 = alpha * beta;
a22 = beta^2;
trA = a11 + a22;
detA = a11*a22 - a12*a12;
disc = max(trA*trA - 4*detA, 0);
emin = (trA - sqrt(disc)) / 2;
emax = (trA + sqrt(disc)) / 2;
end


function scaled = scale_positive_eigen_bounds(bounds, factor_lower, factor_upper)
% Scale positive eigenvalue enclosures by affine perturbation factors.
scaled = I_intval(zeros(size(bounds)));
for k = 1:length(bounds)
    scaled(k) = I_hull(I_inf(factor_lower * bounds(k)), ...
                       I_sup(factor_upper * bounds(k)));
end
end


function val = neumann_complement_coupling(phi_full, psi, P, A_xx, A_yy, A_ux_vy, A_uy_vx)
% Discrete Chat = (P grad phi, J grad psi), equation C-Chat-fem.
% With J grad psi = (-psi_y, psi_x):
%   (P grad phi)_1*(-psi_y) + (P grad phi)_2*psi_x.
val = phi_full' * ( ...
    -P(1,1)*A_ux_vy ...
    -P(1,2)*A_yy ...
    +P(2,1)*A_xx ...
    +P(2,2)*A_uy_vx) * psi;
end


function C = neumann_complement_matrix(P, A_xx, A_yy, A_ux_vy, A_uy_vx)
% Matrix for (P grad phi, J grad psi), used in U-upper-eigenspace.
C = -P(1,1)*A_ux_vy ...
    -P(1,2)*A_yy ...
    +P(2,1)*A_xx ...
    +P(2,2)*A_uy_vx;
end


function cluster_ids = complete_cluster_ids(clusters, cutoff, label)
% Return all cluster ids contained in 1:cutoff; reject partial clusters.
cluster_ids = [];
for c = 1:numel(clusters)
    idx = clusters{c};
    idx = idx(:).';
    if idx(1) > cutoff
        break;
    end
    if idx(end) > cutoff
        error('%s cutoff %d cuts the certified cluster [%d,%d].', ...
            label, cutoff, idx(1), idx(end));
    end
    cluster_ids(end+1) = c; %#ok<AGROW>
end

if isempty(cluster_ids)
    error('%s cutoff %d does not contain any complete cluster.', label, cutoff);
end

last_idx = clusters{cluster_ids(end)};
if last_idx(end) ~= cutoff
    error('%s clusters do not cover exactly 1:%d.', label, cutoff);
end
end


function eta = individual_est_a_from_singleton(lam, lam_h, delta_b_sq)
% Lemma bbar-b-relation for a one-dimensional isolated eigenspace.
db2 = I_sup(delta_b_sq);
if isinf(db2)
    db2 = 1;
end
db2 = min(max(db2, 0), 1);
root_factor = sqrt(max(1 - db2, 0));
eta_sq_upper = I_sup(lam) + I_sup(lam_h) - 2 * I_inf(lam) * root_factor;
eta = I_intval(sqrt(max(eta_sq_upper, 0)));
end


function [proj_lower, proj_upper, coeff] = projection_norm_bounds(reference, basis, coupling_matrix, metric_matrix)
% Certified bounds for sqrt(b^T K^{-1} b) in the revised eigenspace theorem.
%
% The columns are first orthonormalized using the midpoint Gram matrix.  The
% interval Gram error then gives spectral bounds for K and converts the
% coefficient vector norm into lower and upper projection-norm bounds.
[basis_orth, ok, reason] = orthonormalize_wrt_midpoint_metric(basis, metric_matrix);
if ~ok
    error('Cannot orthonormalize cluster basis: %s', reason);
end

coeff = I_intval(zeros(size(basis_orth, 2), 1));
for q = 1:size(basis_orth, 2)
    coeff(q) = reference' * coupling_matrix * basis_orth(:,q);
end

[coeff_norm_lower, coeff_norm_upper] = interval_vector_norm_bounds(coeff);

gram = basis_orth' * metric_matrix * basis_orth;
gram = I_hull(gram, gram');
dim = size(gram, 1);
gram_error = I_sup(norm(I_intval(eye(dim)) - gram, 2));
if ~(gram_error < 1)
    error('Cannot certify positive definiteness of the cluster Gram matrix.');
end

lambda_min_lower = 1 - gram_error;
lambda_max_upper = 1 + gram_error;

proj_lower = I_intval(coeff_norm_lower / sqrt(lambda_max_upper));
proj_upper = I_intval(coeff_norm_upper / sqrt(lambda_min_lower));
end


function val = squared_positive_part_lower(norm_lower, radius)
% Lower bound for (max{norm - radius, 0})^2.
lo = max(I_inf(norm_lower) - I_sup(radius), 0);
val = I_intval(lo^2);
end


function val = squared_sum_upper(norm_upper, radius)
% Upper bound for (norm + radius)^2.
hi = I_sup(norm_upper) + I_sup(radius);
val = I_intval(hi^2);
end


function [energy_lower, proj_h_norm_lower, proj_radius, gram_norm_upper, ok, reason, dhat] = ...
    neumann_cluster_projection_energy_lower(phi_full, psi_block, P, ...
    A_xx, A_yy, A_ux_vy, A_uy_vx, K_N, mu_block, eta_phi_i, eta_psi_cluster, lamih, norm_dotP)
% Certified cluster contribution to the Neumann complement in U.
%
% Paper reference:
%   U_{i,M} = G_i - sum_m |(dotP grad phi_i, J grad psi_m)|^2/mu_m.
%
% For a complete Neumann cluster K, the sum over K is the squared norm of
% the projection of dotP grad phi_i onto
%   J grad span{psi_m : m in K}
% with the basis normalized by the gradient inner product.  This subspace
% norm is invariant under rotations inside a multiple/near-multiple cluster.

energy_lower = I_intval(0);
proj_h_norm_lower = I_intval(0);
proj_radius = I_intval(Inf);
gram_norm_upper = I_intval(Inf);
ok = false;
reason = '';
dhat = I_intval(zeros(size(psi_block, 2), 1));

if isempty(psi_block)
    reason = 'Empty Neumann cluster.';
    return;
end
if isinf(eta_phi_i) || isinf(eta_psi_cluster)
    reason = 'Cluster eigenspace estimator is infinite.';
    return;
end

mu_lower = min(I_inf(mu_block));
if ~(mu_lower > 0)
    reason = 'Cluster Neumann eigenvalue is not provably positive.';
    return;
end

[psi_grad_orth, orth_ok, orth_reason] = orthonormalize_wrt_midpoint_metric(psi_block, K_N);
if ~orth_ok
    reason = orth_reason;
    return;
end

for q = 1:size(psi_grad_orth, 2)
    dhat(q) = neumann_complement_coupling(phi_full, psi_grad_orth(:,q), P, ...
        A_xx, A_yy, A_ux_vy, A_uy_vx);
end

% The columns are orthonormalized with the midpoint matrix.  The interval
% Gram bound converts coefficient-vector norm into a certified projection
% norm on the computed cluster subspace.
gram = psi_grad_orth' * K_N * psi_grad_orth;
gram = I_hull(gram, gram');
gram_norm_upper = I_intval(I_sup(norm(gram, 2)));
if ~(I_sup(gram_norm_upper) > 0) || isinf(I_sup(gram_norm_upper))
    reason = 'Could not certify a finite gradient Gram norm for the cluster basis.';
    return;
end

dhat_norm_lower = interval_vector_norm_lower(dhat);
proj_h_norm_lower = dhat_norm_lower / sqrt(gram_norm_upper);

% Projection error:
%   ||P_K g - P_K^h g_h||
%       <= ||g-g_h|| + gap(J grad E_K, J grad E_K^h) ||g_h||.
% Here g=dotP grad phi_i.  The Liu--Vejchodsky eta is an absolute gradient
% error for the scalar Neumann eigenspace; division by sqrt(mu_min) converts
% it to a gap for the gradient-normalized J grad cluster.
subspace_gap_upper = I_sup(eta_psi_cluster / sqrt(I_intval(mu_lower)));
subspace_gap_upper = min(max(subspace_gap_upper, 0), 1);

phi_error = norm_dotP * eta_phi_i;
subspace_error = I_intval(subspace_gap_upper) * norm_dotP * sqrt(lamih);
proj_radius = phi_error + subspace_error;

lower_norm = max(I_inf(proj_h_norm_lower) - I_sup(proj_radius), 0);
energy_lower = I_intval(lower_norm^2);
ok = true;
end


function [Uo, ok, reason] = orthonormalize_wrt_midpoint_metric(U, M)
% Build a stable basis using the midpoint Gram matrix; interval Gram bounds
% are checked later before the basis is used for certification.
ok = true;
reason = '';

Gint = U' * M * U;
Gsym = (Gint + Gint') / 2;
Gmid = I_mid(Gsym);

[R, pflag] = chol(Gmid);
if pflag ~= 0
    ok = false;
    reason = 'Midpoint gradient Gram matrix is not positive definite.';
    Uo = U;
    return;
end

Uo = U / R;
end


function nrm_lower = interval_vector_norm_lower(v)
% Lower bound for ||v||_2 when each entry is an interval.
[nrm_lower, ~] = interval_vector_norm_bounds(v);
end


function [nrm_lower, nrm_upper] = interval_vector_norm_bounds(v)
% Lower and upper bounds for ||v||_2 when each entry is an interval.
s_lower = 0;
s_upper = 0;
for q = 1:length(v)
    aq = abs(v(q));
    lo = max(I_inf(aq), 0);
    hi = max(I_sup(aq), 0);
    s_lower = s_lower + lo^2;
    s_upper = s_upper + hi^2;
end
nrm_lower = sqrt(s_lower);
nrm_upper = sqrt(s_upper);
end


function [sq_lower, sq_upper] = squared_abs_bounds(center, radius)
% Bounds for (max{|center|-radius,0})^2 and (|center|+radius)^2.
abs_center = abs(center);
lo = max(I_inf(abs_center) - I_sup(radius), 0);
hi = I_sup(abs_center) + I_sup(radius);
sq_lower = I_intval(lo^2);
sq_upper = I_intval(hi^2);
end


function assert_separated(den, label)
if (I_inf(den) <= 0) && (I_sup(den) >= 0)
    error('Eigenvalue enclosures overlap: %s contains 0.', label);
end
end
