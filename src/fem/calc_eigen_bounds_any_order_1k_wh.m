function [eig_bounds, LA_eigf, A_grad, A_L2, A_xx, A_xy, A_yy, ...
    U_with_bdry, meshCG, diagnostics] = ...
    calc_eigen_bounds_any_order_1k_wh( ...
        neig,tri_intval,N_LG,N_rho,LagrangeOrder)

    % =============================================================
    % Step 1: Compute rho <= lambda_{n+1} using CR + explicit C_h
    % =============================================================
    
    mesh_size_rho = 1/N_rho;

    a = tri_intval(5);
    b = tri_intval(6);
        
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
        [lamCR,CR_indices] = ...
            veigs(K_CR,M_CR,neig+1,'sm');
        lamCR = lamCR(:);
        CR_indices = CR_indices(:);
        ordered_CR = I_intval(zeros(neig+1,1));
        for k = 1:neig+1
            location = find(CR_indices == k);
            if numel(location) ~= 1
                error( ...
                    'calc_eigen_bounds_any_order_1k_wh:BadCRIndices', ...
                    ['The verified CR solve did not return precisely ', ...
                     'spectral index %d.'],k);
            end
            ordered_CR(k) = lamCR(location);
        end
        lamCR = ordered_CR;
        Ch = I_intval('0.1893') * hmax;
    else
        lamCR = eigs(K_CR, M_CR, neig + 1, 'sm');
        [~,order_CR] = sort(real(lamCR),'ascend');
        lamCR = lamCR(order_CR);
        Ch = 0.1893 * hmax;
    end

    % rho_j = lam_j / (1 + (Ch^2) * lam_j), then take rho = rho_{n+1}
    rhoCand = lamCR ./ (1 + lamCR .* (Ch^2));
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

    LA_eigf=U; A_grad=K_CG; A_L2=M_CG;

    % The values returned by EIGS above are only floating-point
    % approximations.  Certify the Ritz values of their conforming span
    % before using a number as a min--max upper bound or as the LG shift
    % separation threshold.
    ritzCG = verified_ritz_enclosures(U,K_CG,M_CG,neig);
    Lambda = I_intval(I_sup(ritzCG(neig))); % rigorous upper bound for lambda_n
    if ~(I_sup(Lambda) < I_inf(rho))
        message = ['Separation not verified (need Lambda < rho). ', ...
                   'Refine meshes for CR and/or CG to enforce Lambda < rho.'];
        if INTERVAL_MODE
            error('calc_eigen_bounds_any_order_1k_wh:SeparationFailed', ...
                message);
        else
            warning('calc_eigen_bounds_any_order_1k_wh:SeparationFailed', ...
                message);
        end
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
    if INTERVAL_MODE && ~isspd(B)
        error('calc_eigen_bounds_any_order_1k_wh:BLNotSPD', ...
            'The Lehmann--Goerisch matrix B is not verified positive definite.');
    end

    % Solve A x = mu B x for the smallest neig mu (should be negative)
    mu = I_eig(A, B, neig);

    if any(I_sup(mu) >= 0)
        message = ['Some mu are not verified negative. ', ...
                   'Check separation and assembly.'];
        if INTERVAL_MODE
            error('calc_eigen_bounds_any_order_1k_wh:MuNotNegative', ...
                message);
        else
            warning('calc_eigen_bounds_any_order_1k_wh:MuNotNegative', ...
                message);
        end
    end

    % Lehmann–Goerisch bound:
    %   lambda_i >= rho - rho/(1 - mu_{n+1-i})
    muRev = mu(end:-1:1);
    % The theorem fixes this reverse-index pairing.  Do not reorder
    % overlapping interval enclosures by their midpoints.
    lamLow = (rho-rho./(1-muRev));
    lamLow = lamLow(:);

    eig_bounds = I_intval(zeros(neig,1));
    for k = 1:neig
        lower = I_inf(lamLow(k));
        upper = I_sup(ritzCG(k));
        if lower > upper
            error('calc_eigen_bounds_any_order_1k_wh:InconsistentBounds', ...
                'LG lower bound exceeds the verified Ritz upper bound.');
        end
        eig_bounds(k) = I_infsup(lower,upper);
    end

    % lamCG is deliberately not returned as an upper bound.
    diagnostics = struct();
    diagnostics.method = 'Lehmann-Goerisch-with-CR-Liu-shift';
    diagnostics.CR_discrete_eigenvalues = lamCR;
    diagnostics.CR_Liu_lower_bounds = rhoCand;
    diagnostics.CR_lambda_next_lower = I_intval(I_inf(rho));
    diagnostics.CR_hmax = hmax;
    diagnostics.CR_constant = Ch;
    diagnostics.Ritz_enclosures = ritzCG;
    diagnostics.LG_shift = rho;
    diagnostics.LG_generalized_eigenvalues = mu;
    diagnostics.LG_lower_bounds = lamLow;
    diagnostics.separation_margin = rho-Lambda;
    diagnostics.requires_upper_cluster_dimension = false;

end
