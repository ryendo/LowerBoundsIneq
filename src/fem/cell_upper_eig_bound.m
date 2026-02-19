function cell_ub = cell_upper_eig_bound(cell_data)
% CELL_UPPER_EIG_BOUND: Compute upper bounds for eigenvalues using Lagrange FEM
%
% Paper Reference: Rayleigh-Ritz principle
%   FEM eigenvalues are upper bounds for true eigenvalues
%
% Inputs:
%   cell_data: structure with fields:
%     - x_inf, x_sup, theta_inf, theta_sup: geometry parameters (strings)
%     - fem_order_upper: FEM polynomial order (1 or 2)
%     - mesh_size_upper: mesh size for upper bound computation
%     - neig (optional): number of eigenvalues to compute (default: 1)
%
% Outputs:
%   cell_ub: upper bounds for eigenvalues (vector of length neig)
%
% Author: Based on paper by R. Endo, X. Liu, P. Mariano
% Date: 2025-01-14 (modified for J1/J2 evaluation)

    % Extract and convert geometry parameters to intervals
    x1 = I_intval(str2num(cell_data.x_inf));
    x2 = I_intval(str2num(cell_data.x_sup));
    t1 = I_intval(str2num(cell_data.theta_inf));
    t2 = I_intval(str2num(cell_data.theta_sup));

    % Define vertices
    % Note: a1, b1 corresponds to x1, t1 (use smallest triangle for upper bound)
    a1 = x1;  b1 = x1 * tan(t1);

    % Number of eigenvalues to compute (default: 1 for J evaluation)
    if isfield(cell_data, 'neig')
        neig = cell_data.neig;
    else
        neig = 1;  % Default to lambda_1 only for J1/J2 evaluation
    end

    % Setup FEM parameters
    fem_ord = I_mid(I_intval(cell_data.fem_order_upper));
    mesh_size = cell_data.mesh_size_upper;

    % Mesh Generation (using Gmsh)
    mesh_rho = make_mesh_by_gmsh(a1, b1, mesh_size);
    vert_rho = I_intval(mesh_rho.nodes);
    edge_rho = mesh_rho.edges;
    tri_rho  = mesh_rho.elements;
    bd_rho   = mesh_rho.boundary_edges;

    % Compute upper bounds using Lagrange FEM
    cg_lams = upper_eig_bound(fem_ord, vert_rho, edge_rho, tri_rho, bd_rho, neig);

    % Return the supremum of the requested eigenvalues
    cell_ub = I_sup(cg_lams(1:neig));
end