function [mu_bounds, mu_h, psi_list, K_CG, M_CG, A_xx, A_xy, A_yy, A_ux_vy, A_uy_vx, Est_grad, meshCG] = calc_neumann_eigen_bounds_any_order_1k_wh(neig_positive, tri_intval, N_LG, N_rho, LagrangeOrder, meshCG)
%CALC_NEUMANN_EIGEN_BOUNDS_ANY_ORDER_1K_WH  Neumann data for the DN estimator.
%
% Paper reference:
%   - Neumann weak problem and conforming/CR FEM bounds in Lemma est-tau.
%   - Eigenspace errors eta_{psi,m}^{p,h} used in Theorem main-theorem-est.
%
% neig_positive counts nonconstant Neumann modes:
%   mu_1, ..., mu_neig_positive.
% The zero constant eigenvalue is solved internally and then discarded.

format compact long infsup

if nargin < 6
    meshCG = [];
end

a = I_intval(tri_intval(5));
b = I_intval(tri_intval(6));

% -------------------------------------------------------------------------
% CR lower bounds for positive Neumann eigenvalues, Lemma est-tau.
% The CR Neumann space keeps boundary edge dofs; only the zero mode is later
% discarded from the sorted eigenvalue list.
% -------------------------------------------------------------------------
mesh_size_rho = 1/N_rho;
meshCR = make_mesh_by_gmsh(a, b, mesh_size_rho);

vertCR = meshCR.nodes;
edgeCR = meshCR.edges;
triCR  = meshCR.elements;
triByEdge = find_tri2edge(triCR, edgeCR);

[M_CR, K_CR] = create_matrix_crouzeix_raviart(triCR, edgeCR, vertCR, triByEdge);
hmax = find_mesh_hmax(vertCR, edgeCR);

global INTERVAL_MODE;
if INTERVAL_MODE
    muCR = positive_neumann_cr_veigs(K_CR, M_CR, neig_positive);
    Ch = I_intval('0.1893') * hmax;
else
    muCR = positive_neumann_cr_eigs(K_CR, M_CR, neig_positive);
    Ch = 0.1893 * hmax;
end

muLowCand = muCR ./ (1 + muCR .* (Ch^2));
[~, idx] = sort(I_mid(muLowCand));
muLowCand = muLowCand(idx);

if numel(I_mid(muLowCand)) < neig_positive
    error('CR Neumann eigen computation returned too few eigenvalues.');
end

mu_low = muLowCand(1:neig_positive);

% -------------------------------------------------------------------------
% Conforming Neumann eigenpairs on the same full-dof mesh used by the
% Dirichlet data whenever that mesh is supplied by calc_ddlami_lower_bound.
% -------------------------------------------------------------------------
if isempty(meshCG)
    mesh_size_LG = 1/N_LG;
    meshCG = make_mesh_by_gmsh(a, b, mesh_size_LG);
end

vertCG = meshCG.nodes;
edgeCG = meshCG.edges;
triCG  = meshCG.elements;

[mu_h, psi_list, K_CG, M_CG, A_xx, A_xy, A_yy, A_ux_vy, A_uy_vx] = ...
    laplace_neumann_eig_lagrange_detailed(LagrangeOrder, vertCG, edgeCG, triCG, neig_positive);

mu_h = I_intval(mu_h(:));
mu_low = I_intval(mu_low(:));
mu_bounds = I_hull(mu_h, mu_low);

Est_grad = [];
if nargout >= 11
    clusters = auto_cluster_eigenvalues(mu_bounds, 0.01);
    [Est_grad, ~, ~, ~] = calc_grad_error_bounds(mu_bounds, mu_h, psi_list, K_CG, M_CG, clusters);
end

end


function muCR = positive_neumann_cr_veigs(K_CR, M_CR, neig_positive)
% Certified CR bounds for positive Neumann modes in Lemma est-tau.
% The Neumann CR problem has the exact zero constant mode.  To keep the
% matrices sparse, we do not project it out.  Instead, midpoint eigs provides
% positive spectral shifts, and veigs verifies the CR eigenvalues near those
% shifts.
sigma = positive_neumann_shift_candidates(K_CR, M_CR, neig_positive);

mu_list = I_intval([]);
for s = 1:numel(sigma)
    raw = veigs(K_CR, M_CR, 1, sigma(s));
    raw = raw(:);
    for q = 1:numel(raw)
        if I_sup(raw(q)) > 0
            mu_list(end+1, 1) = raw(q); %#ok<AGROW>
        end
    end
    if numel(I_mid(mu_list)) >= neig_positive
        break;
    end
end

if numel(I_mid(mu_list)) < neig_positive
    error('Could not verify enough positive Neumann CR eigenvalues.');
end

[~, idx] = sort(I_mid(mu_list));
muCR = mu_list(idx(1:neig_positive));
end


function muCR = positive_neumann_cr_eigs(K_CR, M_CR, neig_positive)
% Non-interval fallback: discard the zero Neumann mode from midpoint eigs.
sigma = positive_neumann_shift_candidates(K_CR, M_CR, neig_positive);
muCR = sigma(1:neig_positive);
end


function sigma = positive_neumann_shift_candidates(K_CR, M_CR, neig_positive)
% Positive midpoint CR eigenvalues used only as shifts for veigs.
n = size(K_CR, 1);
request = min(n, max(neig_positive + 4, 8));

while true
    [~, D, flag] = eigs(I_mid(K_CR), I_mid(M_CR), request, 'smallestabs');
    if flag ~= 0
        warning('EIGS did not fully converge for Neumann CR shifts; using available Ritz values.');
    end

    vals = sort(real(diag(D)));
    vals = vals(isfinite(vals));
    scale = max(1, max(abs(vals)));
    zero_tol = 1e-8 * scale;
    sigma = vals(vals > zero_tol);

    if numel(sigma) >= neig_positive
        sigma = sigma(1:neig_positive);
        return;
    end

    if request >= n
        error('Could not find enough positive Neumann CR shift candidates.');
    end
    request = min(n, request + max(neig_positive + 4, 8));
end
end
