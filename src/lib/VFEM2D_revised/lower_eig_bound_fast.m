function eig_bounds = lower_eig_bound_fast(LagrangeOrder, meshCR, meshCG, nEig)
%LOWER_EIG_BOUND_K1  Lehmann–Goerisch lower bounds for the first nEig Dirichlet eigenvalues.
%
% This implements (Algorithm 1) + the construction around (5.25):
%   div w_i^(1) + rho * w_i^(2) + v_i = 0  in Omega,
% with v_i chosen as the FEM eigenfunction u_{i,h}.
%
% INPUTS
%   LagrangeOrder : polynomial degree p for Lagrange FEM (CG space)
%   meshCR        : mesh struct for Crouzeix–Raviart (CR) step (rho bound)
%                  fields: nodes, edges, elements, boundary_edges
%   meshCG        : mesh struct for Lagrange/RT step (Lehmann–Goerisch)
%                  fields: nodes, edges, elements, boundary_edges
%   nEig          : number of eigenvalues to bound (leading nEig)
%
% OUTPUT
%   eig_bounds    : 1-by-nEig vector of guaranteed lower bounds (if INTERVAL_MODE==true,
%                  returns the lower endpoints of interval bounds)

    % =============================================================
    % Step 1: Compute rho <= lambda_{n+1} using CR + explicit C_h
    % =============================================================
    vertCR = meshCR.nodes;
    edgeCR = meshCR.edges;
    triCR  = meshCR.elements;
    bdCR   = meshCR.boundary_edges;

    isBnd = ismember(edgeCR, bdCR, 'rows');
    triByEdge = find_tri2edge(triCR, edgeCR);
    bdEdgeIds = find(isBnd > 0);

    [M_CR, K_CR] = create_matrix_crouzeix_raviart(triCR, edgeCR, vertCR, triByEdge);

    nEdges = size(edgeCR, 1);
    dof = 1:nEdges;
    dof(bdEdgeIds) = [];

    M_CR = M_CR(dof, dof);
    K_CR = K_CR(dof, dof);

    hmax = find_mesh_hmax(vertCR, edgeCR);

    % Solve CR generalized eigenproblem: K u = lambda M u
    global INTERVAL_MODE;
    if INTERVAL_MODE
        lamCR = veigs(K_CR, M_CR, nEig + 1, 'sm');
        Ch = I_intval('0.1893') * hmax;
    else
        lamCR = eigs(K_CR, M_CR, nEig + 1, 'sm');
        Ch = 0.1893 * hmax;
    end

    % rho_j = lam_j / (1 + (Ch^2) * lam_j), then take rho = rho_{n+1}
    rhoCand = lamCR ./ (1 + lamCR .* (Ch^2));
    [~, idx] = sort(I_mid(rhoCand));
    rhoCand = rhoCand(idx);

    nCand = numel(I_mid(rhoCand));
    if nCand < nEig + 1
        error('CR eigen computation returned too few eigenvalues.');
    end
    rho = rhoCand(nEig + 1);

    % =============================================================
    % Step 2: Compute CG (Lagrange) eigenpairs (upper bounds)
    % =============================================================
    vertCG = meshCG.nodes;
    edgeCG = meshCG.edges;
    triCG  = meshCG.elements;
    bdCG   = meshCG.boundary_edges;

    % Expected outputs (as in your original code):
    %   lamCG: nEig-by-1 approximate eigenvalues (upper bounds)
    %   U    : (dofs-by-nEig) eigenfunctions on interior dofs
    %   Ubd  : (full-dofs-by-nEig) eigenfunctions including boundary dofs
    %   K,M  : assembled stiffness/mass on the same interior dof basis as U
    [lamCG, U, U_with_bdry, K_CG, M_CG, ~, ~, ~, ~] = ...
        laplace_eig_lagrange_detailed(LagrangeOrder, vertCG, edgeCG, triCG, bdCG, nEig);

    Lambda = max(lamCG(1:nEig)); % upper bound for lambda_n
    if ~(I_sup(Lambda) < I_inf(rho))
        warning(['Separation not verified (need Lambda < rho). ', ...
                 'Refine meshes for CR and/or CG to enforce Lambda < rho.']);
    end

    % =============================================================
    % Step 3: Build w_i = {w_i^(1), w_i^(2)}
    % =============================================================
    % IMPORTANT:
    % RT_Hdiv_problem_dirichlet returns:
    %   W   : coefficient matrix whose i-th column represents w_i in X
    %   BX  : matrix representing b_G(.,.) on X so that A2 = W' * BX * W
    %    
    % lambda_hat = I_mid(rho);
    lambda_hat = 1e-14;
    rhoLG = rho + lambda_hat;

    RTorder = LagrangeOrder;
    [W, BX] = RT_Hdiv_problem_dirichlet_fast( ...
        RTorder, LagrangeOrder, vertCG, edgeCG, triCG, bdCG, U_with_bdry, lambda_hat);

    % =============================================================
    % Step 4: Assemble LG matrices and solve interval generalized EVP
    % =============================================================    

    A0 = U' * (K_CG + lambda_hat * M_CG) * U;
    A1 = U' * M_CG * U;
    A2 = W' * BX * W;     % b_G(w_i, w_j)

    A = A0 - rhoLG * A1;
    B = A0 - 2*rhoLG * A1 + (rhoLG*rhoLG) * A2;

    % Symmetrize (interval I_hull if interval mode)
    A = I_hull(A, A');
    B = I_hull(B, B');

    % Solve A x = mu B x for the smallest nEig mu (should be negative)
    if INTERVAL_MODE
        mu = I_eig(A, B, nEig);
    else
        muAll = eig(A, B);
        muAll = sort(real(muAll), 'ascend');
        mu = muAll(1:nEig);
    end
    mu = mu(:);

    if any(I_sup(mu) >= 0)
        warning('Some mu are not verified negative. Check separation and assembly.');
    end

    % Lehmann–Goerisch bound:
    %   lambda_i >= rho - rho/(1 - mu_{n+1-i})
    muRev = mu(end:-1:1);
    lamLow_shift = rhoLG - rhoLG ./ (1 - muRev);
    lamLow = lamLow_shift - lambda_hat;

    % Sort by I_midpoints and return verified lower endpoints
    [~, idx] = sort(I_mid(lamLow));
    lamLow = lamLow(idx);

    eig_bounds = (I_inf(lamLow))'; % 1-by-nEig

end