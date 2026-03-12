function [eig_bounds, U, K_CG, M_CG, A_xx, A_xy, A_yy] = calc_eigen_bounds_any_order_1k_wh(neig,tri_intval,N_LG,N_rho,LagrangeOrder)

    % =============================================================
    % Step 1: Compute rho <= lambda_{n+1} using CR + explicit C_h
    % =============================================================
    format compact long infsup
    mesh_size_rho = 1/N_rho;

    a = I_intval(tri_intval(5));
    b = I_intval(tri_intval(6));
        
    meshCR = make_mesh_by_gmsh(a, b, mesh_size_rho);

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
        lamCR = veigs(K_CR, M_CR, neig + 1, 'sm');
        Ch = I_intval('0.1893') * hmax;
    else
        lamCR = eigs(K_CR, M_CR, neig + 1, 'sm');
        Ch = 0.1893 * hmax;
    end

    % rho_j = lam_j / (1 + (Ch^2) * lam_j), then take rho = rho_{n+1}
    rhoCand = lamCR ./ (1 + lamCR .* (Ch^2));
    [~, idx] = sort(I_mid(rhoCand));
    rhoCand = rhoCand(idx);
    
    nCand = numel(I_mid(rhoCand));
    if nCand < neig + 1
        error('CR eigen computation returned too few eigenvalues.');
    end
    rho = rhoCand(neig + 1);

    % =============================================================
    % Step 2: Compute CG (Lagrange) eigenpairs (upper bounds)
    % =============================================================

    mesh_size_LG = 1/N_LG;
    meshCG = make_mesh_by_gmsh(a, b, mesh_size_LG);

    vertCG = meshCG.nodes;
    edgeCG = meshCG.edges;
    triCG  = meshCG.elements;
    bdCG   = meshCG.boundary_edges;

    % Expected outputs (as in your original code):
    %   lamCG: neig-by-1 approximate eigenvalues (upper bounds)
    %   U    : (dofs-by-neig) eigenfunctions on interior dofs
    %   Ubd  : (full-dofs-by-neig) eigenfunctions including boundary dofs
    %   K,M  : assembled stiffness/mass on the same interior dof basis as U

    [lamCG, U, U_with_bdry, K_CG, M_CG, A_xx, A_xy, A_yy, ~] = ...    
        laplace_eig_lagrange_detailed(LagrangeOrder, vertCG, edgeCG, triCG, bdCG, neig);

    Lambda = max(lamCG(1:neig)); % upper bound for lambda_n
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
    RTorder = LagrangeOrder;
    A2 = RT_Hdiv_problem_dirichlet(meshCG,RTorder,U_with_bdry);

    % =============================================================
    % Step 4: Assemble LG matrices and solve interval generalized EVP
    % =============================================================
    A0 = U' * K_CG * U;   % (grad v_i, grad v_j)
    A1 = U' * M_CG * U;   % (v_i, v_j)

    A = A0 - rho * A1;
    B = A0 - 2*rho * A1 + (rho*rho) * A2;

    % Symmetrize (interval I_hull if interval mode)
    A = I_hull(A, A');
    B = I_hull(B, B');

    % Solve A x = mu B x for the smallest neig mu (should be negative)
    mu = I_eig(A, B, neig);

    if any(I_sup(mu) >= 0)
        warning('Some mu are not verified negative. Check separation and assembly.');
    end

    % Lehmann–Goerisch bound:
    %   lambda_i >= rho - rho/(1 - mu_{n+1-i})
    muRev = mu(end:-1:1);
    lamLow = rho - rho ./ (1 - muRev);

    % Sort by I_midpoints and return verified lower endpoints
    [~, idx] = sort(I_mid(lamLow));
    lamLow = I_intval(I_inf(lamLow(idx)'));

    eig_bounds = I_hull(lamCG,lamLow); % 1-by-neig

    % lamCG
    % lamLow
    % eig_bounds

end