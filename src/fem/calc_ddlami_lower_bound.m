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
    M_neumann = N_spectral;
end

if i > N_spectral
    error('The paper bound requires i <= N_spectral.');
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

x_base   = base_triangle(5);
y_base   = base_triangle(6);
x_target = triangle(5);
y_target = triangle(6);

% Eigenvalue perturbation over the target cell, used exactly as in the
% paper's uniform cell algorithms for Omega_up.
S_inv     = [1, (x_base - x_target) / y_target; 0, y_base / y_target];
SinvSinvt = S_inv * S_inv';
scale_interval = norm(SinvSinvt);

lams_interval = lams_raw   * scale_interval;
lams_h        = lams_h_raw * scale_interval;

clusters = auto_cluster_eigenvalues(lams_interval, 0.01);
eta_phi = calc_grad_error_bounds(lams_interval, lams_h, uh_list, K_D, M_D, clusters);

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
mu_interval = mu_raw * scale_interval;
mu_h = mu_h_raw * scale_interval;
mu_clusters = auto_cluster_eigenvalues(mu_interval, 0.01);
[eta_psi, eta_psi_cluster_sq, ~, eta_psi_info] = ...
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
e_dlami = norm_dotP * (sqrt(lami) + sqrt(lamih)) * eta_phi(i);
dlami = I_hull(dlamih - e_dlami, dlamih + e_dlami);

% -------------------------------------------------------------------------
% A_i and epsilon_A: equations explicit-eps-A and explicit-A-bounds.
% -------------------------------------------------------------------------
Ahat_i = uh_list(:,i)' * shape_D2 * uh_list(:,i);
eps_A = norm_ddotP * eta_phi(i) * (sqrt(lami) + sqrt(lamih));
A_lower = Ahat_i - eps_A;
A_upper = Ahat_i + eps_A;

% -------------------------------------------------------------------------
% B_ik, epsilon_B, underline/overline mathcal B: equations
% explicit-eps-B and explicit-Bsq-bounds.
% -------------------------------------------------------------------------
Bhat = I_intval(zeros(N_spectral, 1));
Bsq_lower = I_intval(zeros(N_spectral, 1));
Bsq_upper = I_intval(zeros(N_spectral, 1));
for k = 1:N_spectral
    Bhat(k) = uh_list(:,i)' * shape_D * uh_list(:,k);
    eps_B = norm_dotP * (eta_phi(i) * sqrt(lams_interval(k)) + sqrt(lamih) * eta_phi(k));
    [Bsq_lower(k), Bsq_upper(k)] = squared_abs_bounds(Bhat(k), eps_B);
end

% -------------------------------------------------------------------------
% underline D_{i,N}^{p,h}: equation explicit-D-lower-def.
% The sign of lambda_k-lambda_i determines whether the lower or upper square
% bound is used.
% -------------------------------------------------------------------------
D_lower = A_lower;
D_upper = A_upper;
for k = 1:N_spectral
    if k == i
        continue;
    end

    den = lams_interval(k) - lams_interval(i);
    assert_separated(den, sprintf('lambda_%d - lambda_%d', k, i));

    if k < i
        D_lower = D_lower - 2 * Bsq_lower(k) / den;
        D_upper = D_upper - 2 * Bsq_upper(k) / den;
    else
        D_lower = D_lower - 2 * Bsq_upper(k) / den;
        D_upper = D_upper - 2 * Bsq_lower(k) / den;
    end
end

% -------------------------------------------------------------------------
% underline S_{i,N}^{p,h}: equation explicit-S-lower-def.  Unlike D, this
% sum includes k=i; for N=1 this is the first-derivative contribution.
% -------------------------------------------------------------------------
S_lower = I_intval(0);
for k = 1:N_spectral
    if I_inf(lams_interval(k)) <= 0
        error('lambda_%d is not provably positive.', k);
    end
    S_lower = S_lower + Bsq_lower(k) / lams_interval(k);
end

% -------------------------------------------------------------------------
% overline U_{i,M}^{p,h}: equation explicit-U-upper-def.
% G is ||dotP grad phi_i||^2.  The paper writes the Neumann complement as a
% mode-by-mode sum of |C_im|^2/mu_m.  In clustered eigenvalues the individual
% eigenfunctions may rotate, so we certify the invariant cluster quantity
%
%   || P_{J grad E_K^N} (dotP grad phi_i) ||^2
%     = sum_{m in K} |C_im|^2/mu_m,
%
% for every complete Neumann cluster K contained in {1,...,M}.  This is the
% same term in U, but evaluated as a subspace projection instead of as
% unstable individual coefficients.
% -------------------------------------------------------------------------
Ghat = gradient_image_norm_sq(uh_list(:,i), dotP, A_xx, A_xy, A_yy);
eps_G = norm_dotP^2 * eta_phi(i) * (sqrt(lami) + sqrt(lamih));
G_upper = Ghat + eps_G;

U_upper = G_upper;
C_hat = I_intval(zeros(M_neumann, 1));
eps_C_list = I_intval(zeros(M_neumann, 1));
Csq_lower_list = I_intval(zeros(M_neumann, 1));
for m = 1:M_neumann
    Chat = neumann_complement_coupling(uh_list_full(:,i), psi_list(:,m), dotP, N_xx, N_yy, N_ux_vy, N_uy_vx);
    eps_C = norm_dotP * (eta_phi(i) * sqrt(mu_interval(m)) + sqrt(lamih) * eta_psi(m));
    [Csq_lower, ~] = squared_abs_bounds(Chat, eps_C);
    C_hat(m) = Chat;
    eps_C_list(m) = eps_C;
    Csq_lower_list(m) = Csq_lower;
end

num_mu_clusters = numel(mu_clusters);
cluster_energy_lower = I_intval(zeros(num_mu_clusters, 1));
cluster_proj_h_norm_lower = I_intval(zeros(num_mu_clusters, 1));
cluster_proj_radius = I_intval(Inf(num_mu_clusters, 1));
cluster_gram_norm_upper = I_intval(Inf(num_mu_clusters, 1));
cluster_used = false(num_mu_clusters, 1);
cluster_reason = repmat({''}, num_mu_clusters, 1);
cluster_hat_coeff = cell(num_mu_clusters, 1);

for c = 1:num_mu_clusters
    idx_cluster = mu_clusters{c};
    idx_cluster = idx_cluster(:).';

    if idx_cluster(1) > M_neumann
        cluster_reason{c} = 'Cluster starts after the requested Neumann cutoff M.';
        continue;
    end
    if idx_cluster(end) > M_neumann
        cluster_reason{c} = 'Cluster intersects the cutoff M; skipped to avoid subtracting a partial eigenspace.';
        continue;
    end

    if any(I_inf(mu_interval(idx_cluster)) <= 0)
        error('A Neumann cluster used in U contains a non-positive eigenvalue lower bound.');
    end

    eta_cluster = sqrt(eta_psi_cluster_sq(c));
    [Ek_lower, proj_h_norm_lower, proj_radius, gram_norm_upper, ok, reason, dhat] = ...
        neumann_cluster_projection_energy_lower( ...
            uh_list_full(:,i), psi_list(:,idx_cluster), dotP, ...
            N_xx, N_yy, N_ux_vy, N_uy_vx, K_N, ...
            mu_interval(idx_cluster), eta_phi(i), eta_cluster, lamih, norm_dotP);

    cluster_energy_lower(c) = Ek_lower;
    cluster_proj_h_norm_lower(c) = proj_h_norm_lower;
    cluster_proj_radius(c) = proj_radius;
    cluster_gram_norm_upper(c) = gram_norm_upper;
    cluster_used(c) = ok;
    cluster_reason{c} = reason;
    cluster_hat_coeff{c} = dhat;

    if ok
        U_upper = U_upper - Ek_lower;
    end
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
diagnostics.eta_phi_i = eta_phi(i);
diagnostics.eta_psi = eta_psi(1:M_neumann);
diagnostics.eta_psi_cluster = sqrt(eta_psi_cluster_sq);
diagnostics.eta_psi_cluster_sq = eta_psi_cluster_sq;
diagnostics.eta_psi_info = eta_psi_info;
diagnostics.mu_clusters = mu_clusters;
diagnostics.mu_interval = mu_interval(1:M_neumann);
diagnostics.mu_h = mu_h(1:M_neumann);
diagnostics.C_hat = C_hat;
diagnostics.eps_C = eps_C_list;
diagnostics.Csq_lower = Csq_lower_list;
diagnostics.cluster_energy_lower = cluster_energy_lower;
diagnostics.cluster_proj_h_norm_lower = cluster_proj_h_norm_lower;
diagnostics.cluster_proj_radius = cluster_proj_radius;
diagnostics.cluster_gram_norm_upper = cluster_gram_norm_upper;
diagnostics.cluster_used = cluster_used;
diagnostics.cluster_reason = cluster_reason;
diagnostics.cluster_hat_coeff = cluster_hat_coeff;

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
Gmid = mid(Gsym);

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
s = 0;
for q = 1:numel(v)
    aq = abs(v(q));
    lo = max(I_inf(aq), 0);
    s = s + lo^2;
end
nrm_lower = I_intval(sqrt(s));
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
