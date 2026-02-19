function [eig_bounds, LA_eigf, A_grad, A_L2, A_xx, A_xy, A_yy] = calc_eigen_bounds_any_order_1k_wh(neig,tri_intval,N_LG,N_rho,ord,ord_LG)

    mesh_size_LG = 1/N_LG;
    mesh_size_rho = 1/N_rho;
    
    % Step 1: Preliminary CR calculation to get shift parameter 'rho'
    % -------------------------------------------------------------
    
    % Geometry for mesh generation
    a = tri_intval(5);
    b = tri_intval(6);
    
    mesh_rho = make_mesh_by_gmsh(a, b, mesh_size_rho);
    vert_rho = mesh_rho.nodes;
    edge_rho = mesh_rho.edges;
    tri_rho  = mesh_rho.elements;
    bd_rho   = mesh_rho.boundary_edges;
    is_bnd = ismember(edge_rho, bd_rho, 'rows');

    global INTERVAL_MODE
    tri_by_edge = find_tri2edge(tri_rho, edge_rho);
    bd_edge_ids = find(is_bnd > 0);
    
    [A0, A1] = create_matrix_crouzeix_raviart(tri_rho, edge_rho, vert_rho, tri_by_edge);
    
    ne = size(edge_rho, 1);
    dof_idx = 1:ne;
    dof_idx(bd_edge_ids) = [];
    
    CR_A0 = A0(dof_idx, dof_idx);
    CR_A1 = A1(dof_idx, dof_idx);
    
    hmax = find_mesh_hmax(vert_rho, edge_rho);
    

    if INTERVAL_MODE
        CR_eig = veigs(CR_A1, CR_A0, neig, 'sm');
    else
        CR_eig = eigs(CR_A1, CR_A0, neig, 'sm');
    end
    
    Ch_val = I_intval('0.1893') * hmax;
    
    % CR-based lower bounds
    temp_bounds = CR_eig ./ (1 + CR_eig .* (Ch_val^2));
    [~, idx] = sort(I_mid(temp_bounds));
    temp_bounds = temp_bounds(idx);

    rho = temp_bounds(neig);
    
    % Step 2: High-order Lehmann-Goerisch Setup
    % -------------------------------------------------------------
    % Define Triangle for LG
    tri_intval = [I_intval('0'), I_intval('0'); ...
                    I_intval('1'), I_intval('0'); ...
                    a,             b];
                    
    % Make mesh for LG
    mesh_LG = make_mesh_by_gmsh(a, b, mesh_size_LG);
    vert_LG = mesh_LG.nodes;
    edge_LG = mesh_LG.edges;
    tri_LG  = mesh_LG.elements;
    bd_LG   = mesh_LG.boundary_edges;
    is_bnd_lg = ismember(edge_LG, bd_LG, 'rows');
    
    
    Lagrange_order = ord_LG;

    % Compute High-Order Upper Bounds & Eigenfunctions (Lagrange)
    [LA_eig, LA_eigf, LA_eigf_with_bdry, LA_A, LA_M, A_xx, A_xy, A_yy, ~] = ...
        laplace_eig_lagrange_detailed(Lagrange_order, vert_LG, edge_LG, tri_LG, bd_LG, neig);

    A_grad = LA_A; A_L2 = LA_M;

    % Step 3: Construct LG Matrices
    % -------------------------------------------------------------
    A0 = LA_eigf' * LA_A * LA_eigf;
    A1 = LA_eigf' * LA_M * LA_eigf;
    
    RT_order = Lagrange_order;
    
    % Solve H(div) problem using Raviart-Thomas (RT) elements
    [mat_pih, RTRT] = RT_Hdiv_problem_dirichlet(RT_order, Lagrange_order, vert_LG, edge_LG, tri_LG, bd_LG, LA_eigf_with_bdry);

    A2 = mat_pih' * RTRT * mat_pih;
    
    % Matrices for Generalized Eigenvalue Problem
    AL = A0 - rho * A1;
    BL = A0 - 2 * rho * A1 + rho * rho * A2;
    
    % Step 4: Solve Interval Generalized Eigenvalue Problem
    % -------------------------------------------------------------
    % Symmetrize via Interval Hull
    AL = I_hull(AL, AL'); 
    BL = I_hull(BL, BL');

    mus = I_eig(AL, BL, neig);

    % Apply Lehmann-Goerisch Theorem
    LG_eig_low = rho - rho ./ (1 - mus(end:-1:1));
    
    % Sort and extract
    [~, idx] = sort(I_mid(LG_eig_low));
    LG_eig_low = LG_eig_low(idx);
    
    eig_bounds_ = I_intval(LG_eig_low);
    [~, idx] = sort(I_mid(eig_bounds_));
    eig_bounds = eig_bounds_(idx)';
 
end

