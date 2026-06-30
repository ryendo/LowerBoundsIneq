function [eig_value, eig_func, A_grad, A_L2, A_xx, A_xy, A_yy, A_ux_vy, A_uy_vx] = laplace_neumann_eig_lagrange_detailed(lagrange_order, vert, edge, tri, neig_positive)
%LAPLACE_NEUMANN_EIG_LAGRANGE_DETAILED  Conforming Neumann FEM eigenpairs.
%
% Paper reference:
%   Section "Finite element approximation of Dirichlet and Neumann data",
%   equations (Neumann weak problem) and (cg-problem) for V_h^N.
%
% The returned eigenpairs are the first neig_positive nonconstant Neumann
% modes.  The zero constant mode is computed only so it can be discarded.
% Matrices are assembled on the full H^1 Lagrange space, so eig_func columns
% are compatible with Dirichlet eigenfunctions embedded with boundary dofs.
%
% The nonsymmetric derivative matrices A_ux_vy and A_uy_vx are returned for
% the paper quantity
%   C_im = (dotP grad phi_i, J grad psi_m),
% where J = [0 -1; 1 0].

ne = size(edge, 1);
nt = size(tri,  1);
[tri, tri2edge] = local_find_tri2edge(edge, tri, ne, nt);

[A_grad, A_L2, A_xx, A_xy, A_yy, A_ux_vy, A_uy_vx] = ...
    Lagrange_full_eig_matrix(lagrange_order, vert, edge, tri, tri2edge);

ndof = size(A_grad, 1);
if neig_positive >= ndof
    error('Requested too many nonconstant Neumann eigenpairs for this mesh.');
end

request = min(ndof, max(neig_positive + 1, 2));
while true
    [V, D] = eigs(I_mid(A_grad), I_mid(A_L2), request, 'sm');
    vals = real(diag(D));
    [vals, idx] = sort(vals);
    V = real(V(:, idx));

    scale = max(1, max(abs(vals)));
    zero_tol = 1e-8 * scale;
    pos_idx = find(vals > zero_tol);

    if numel(pos_idx) >= neig_positive
        pos_idx = pos_idx(1:neig_positive);
        eig_value = vals(pos_idx);
        eig_func = V(:, pos_idx);
        break;
    end

    if request >= ndof
        error('Could not extract enough nonconstant Neumann eigenpairs.');
    end
    request = min(ndof, request + neig_positive + 1);
end

% Choose the mean-zero representative modulo constants and normalize in the
% L2 inner product used in the Neumann Lehmann--Goerisch theorem.
Mmid = I_mid(A_L2);
one_vec = ones(size(eig_func, 1), 1);
area_mass = one_vec' * Mmid * one_vec;
for k = 1:size(eig_func, 2)
    mean_coeff = (one_vec' * Mmid * eig_func(:,k)) / area_mass;
    eig_func(:,k) = eig_func(:,k) - mean_coeff * one_vec;
    nrm = sqrt(eig_func(:,k)' * Mmid * eig_func(:,k));
    eig_func(:,k) = eig_func(:,k) / nrm;
end

end


function [A_glob_grad, A_glob_L2, A_glob_ux_ux, A_glob_ux_uy_sym, A_glob_uy_uy, A_glob_ux_vy, A_glob_uy_vx] = Lagrange_full_eig_matrix(lagrange_order, vert, edge, tri, tri2edge)
[basis, nbasis] = Lagrange_basis(lagrange_order);
M_ip_elem = Lagrange_inner_product_L1L2L3_all(lagrange_order);

grad_basis = cell(nbasis, 1);
scalar_basis = cell(nbasis, 1);
for i = 1:nbasis
    grad_basis{i} = Lagrange_create_coord_basis_grad(basis, i, lagrange_order);
    scalar_basis{i} = Lagrange_create_coord_basis(basis, i, lagrange_order);
end

M_ip_basis_ij = I_zeros(nbasis, nbasis);
M_ip_basis_grad_ijT_all = cell(nbasis, nbasis);
for i = 1:nbasis
    for j = 1:nbasis
        M_ip_basis_grad_ijT_all{i, j} = grad_basis{i}' * M_ip_elem * grad_basis{j};
        M_ip_basis_ij(i, j) = scalar_basis{i}' * M_ip_elem * scalar_basis{j};
    end
end

nv = size(vert, 1);
ne = size(edge, 1);
nt = size(tri,  1);

ndof = nv + (lagrange_order-1)*ne + (lagrange_order-1)*(lagrange_order-2)/2*nt;

A_glob_grad     = I_sparse(ndof, ndof);
A_glob_L2       = I_sparse(ndof, ndof);
A_glob_ux_ux    = I_sparse(ndof, ndof);
A_glob_ux_vy    = I_sparse(ndof, ndof);
A_glob_uy_vx    = I_sparse(ndof, ndof);
A_glob_uy_uy    = I_sparse(ndof, ndof);

for k = 1:nt
    vert_idx = tri(k,:);
    x1 = vert(vert_idx(1), 1); y1 = vert(vert_idx(1), 2);
    x2 = vert(vert_idx(2), 1); y2 = vert(vert_idx(2), 2);
    x3 = vert(vert_idx(3), 1); y3 = vert(vert_idx(3), 2);

    B = [x2-x1, x3-x1; y2-y1, y3-y1];
    Binv = [y3-y1, x1-x3; y1-y2, x2-x1] / det(B);

    A_grad  = I_zeros(nbasis, nbasis);
    A_ux_vx = I_zeros(nbasis, nbasis);
    A_ux_vy = I_zeros(nbasis, nbasis);
    A_uy_vx = I_zeros(nbasis, nbasis);
    A_uy_vy = I_zeros(nbasis, nbasis);

    for i = 1:nbasis
        for j = 1:nbasis
            mat_base = Binv' * M_ip_basis_grad_ijT_all{i, j} * Binv;
            A_grad(i, j)  = trace(mat_base) * det(B);
            A_ux_vx(i, j) = mat_base(1,1) * det(B);
            A_ux_vy(i, j) = mat_base(1,2) * det(B);
            A_uy_vx(i, j) = mat_base(2,1) * det(B);
            A_uy_vy(i, j) = mat_base(2,2) * det(B);
        end
    end

    A_L2 = M_ip_basis_ij * det(B);

    edge_idx = tri2edge(k, :);
    local_edge_start_vert = [2, 3, 1];
    map_dof_idx_l2g = zeros(nbasis, 1);
    map_dof_idx_l2g(1:3) = vert_idx;
    for i = 1:3
        if edge(edge_idx(i), 1) == vert_idx(local_edge_start_vert(i))
            map_dof_idx_l2g(3+(lagrange_order-1)*(i-1)+1:3+(lagrange_order-1)*i) = ...
                nv+(lagrange_order-1)*(edge_idx(i)-1)+1:1:nv+(lagrange_order-1)*edge_idx(i);
        else
            map_dof_idx_l2g(3+(lagrange_order-1)*(i-1)+1:3+(lagrange_order-1)*i) = ...
                nv+(lagrange_order-1)*edge_idx(i):-1:nv+(lagrange_order-1)*(edge_idx(i)-1)+1;
        end
    end
    map_dof_idx_l2g(3+(lagrange_order-1)*3+1:end) = nv + (lagrange_order-1)*ne + ...
        ((lagrange_order-1)*(lagrange_order-2)/2*(k-1)+1:(lagrange_order-1)*(lagrange_order-2)/2*k);

    A_glob_grad(map_dof_idx_l2g, map_dof_idx_l2g)  = A_glob_grad(map_dof_idx_l2g, map_dof_idx_l2g)  + A_grad;
    A_glob_ux_ux(map_dof_idx_l2g, map_dof_idx_l2g) = A_glob_ux_ux(map_dof_idx_l2g, map_dof_idx_l2g) + A_ux_vx;
    A_glob_ux_vy(map_dof_idx_l2g, map_dof_idx_l2g) = A_glob_ux_vy(map_dof_idx_l2g, map_dof_idx_l2g) + A_ux_vy;
    A_glob_uy_vx(map_dof_idx_l2g, map_dof_idx_l2g) = A_glob_uy_vx(map_dof_idx_l2g, map_dof_idx_l2g) + A_uy_vx;
    A_glob_uy_uy(map_dof_idx_l2g, map_dof_idx_l2g) = A_glob_uy_uy(map_dof_idx_l2g, map_dof_idx_l2g) + A_uy_vy;
    A_glob_L2(map_dof_idx_l2g, map_dof_idx_l2g)    = A_glob_L2(map_dof_idx_l2g, map_dof_idx_l2g)    + A_L2;
end

% Existing Dirichlet code uses A_xy = 1/2 int(u_x v_y + u_y v_x).
A_glob_ux_uy_sym = (A_glob_ux_vy + A_glob_uy_vx) / 2;
end


function M_ip = Lagrange_inner_product_L1L2L3_all(lagrange_order)
ijk = create_ijk(lagrange_order);
len = size(ijk, 1);
M_ip = I_zeros(len, len);
for p = 1:len
    for q = p:len
        pi = ijk(p, 1);
        pj = ijk(p, 2);
        pk = ijk(p, 3);
        qi = ijk(q, 1);
        qj = ijk(q, 2);
        qk = ijk(q, 3);
        M_ip(p, q) = Lagrange_integral_L1L2L3_ijk(pi+qi, pj+qj, pk+qk);
    end
end
M_ip = M_ip + triu(M_ip, 1)';
end


function e = Lagrange_create_coord_basis(basis, idx, Lgrange_order)
len = Lagrange_get_nbasis(Lgrange_order);
e = zeros(len, 1);

i = basis(idx, 1);
j = basis(idx, 2);
k = basis(idx, 3);
e(map_ijk_to_idx(i, j, k, Lgrange_order), 1) = 1;
end


function e = Lagrange_create_coord_basis_grad(basis, idx, Lgrange_order)
len = Lagrange_get_nbasis(Lgrange_order);
dudx = zeros(len, 1);
dudy = zeros(len, 1);

i = basis(idx, 1);
j = basis(idx, 2);
k = basis(idx, 3);

dudx(map_ijk_to_idx(i,   j,   k,   Lgrange_order), 1) = -i + j;
dudx(map_ijk_to_idx(i-1, j+1, k,   Lgrange_order), 1) = -i;
dudx(map_ijk_to_idx(i-1, j,   k+1, Lgrange_order), 1) = -i;
dudx(map_ijk_to_idx(i+1, j-1, k,   Lgrange_order), 1) =  j;
dudx(map_ijk_to_idx(i,   j-1, k+1, Lgrange_order), 1) =  j;

dudy(map_ijk_to_idx(i,   j,   k,   Lgrange_order), 1) = -i + k;
dudy(map_ijk_to_idx(i-1, j+1, k,   Lgrange_order), 1) = -i;
dudy(map_ijk_to_idx(i-1, j,   k+1, Lgrange_order), 1) = -i;
dudy(map_ijk_to_idx(i+1, j,   k-1, Lgrange_order), 1) =  k;
dudy(map_ijk_to_idx(i,   j+1, k-1, Lgrange_order), 1) =  k;

e = [dudx, dudy];
end


function idx = map_ijk_to_idx(i, j, k, n)
idx = (n-i)*(n-i+1)/2 + (n-i-j) + 1;
if i + j + k ~= n
    idx = -1;
end
if i<0 || j<0 || k<0
    idx = [];
end
end


function ijk = create_ijk(n)
ijk = zeros((n+1)*(n+2)/2, 3);
index = 1;
for p = n : -1 : 0
    for q = n-p : -1 : 0
        ijk(index, :) = [p, q, n-p-q];
        index = index + 1;
    end
end
end


function [basis, nbasis] = Lagrange_basis(lagrange_order)
n      = lagrange_order;
nbasis = Lagrange_get_nbasis(lagrange_order);
basis  = zeros(nbasis, 3);

index = 1;
basis(index, :) = [n, 0, 0];
index = index + 1;
if n > 0
    basis(index, :) = [0, n, 0];
    index = index + 1;
    basis(index, :) = [0, 0, n];
    index = index + 1;
end
for p = n-1 : -1 : 1
    basis(index, :) = [0, p, n-p];
    index = index + 1;
end
for p = n-1 : -1 : 1
    basis(index, :) = [n-p, 0, p];
    index = index + 1;
end
for p = n-1 : -1 : 1
    basis(index, :) = [p, n-p, 0];
    index = index + 1;
end
for p = n-2 : -1 : 1
    for q = n-1-p : -1 : 1
        basis(index, :) = [p, q, n-p-q];
        index = index + 1;
    end
end
end


function nbasis = Lagrange_get_nbasis(lagrange_order)
nbasis = (lagrange_order+1) * (lagrange_order+2) / 2;
end


function y = Lagrange_integral_L1L2L3_ijk(i, j, k)
y = I_intval(factorial(i) * factorial(j) * factorial(k)) / I_intval(factorial(i+j+k+2));
end


function [tri, tri2edge] = local_find_tri2edge(edge, tri, ne, nt)
tri2edge = zeros(nt, 3);
edge_idx = sort(edge, 2) * [ne; 1];
for k = 1:nt
    edge_local = [2 3; 1 3; 1 2];
    value = sort(reshape(tri(k, edge_local), 3, 2), 2) *[ne; 1];
    [~, idx] = ismember(value, edge_idx);
    tri2edge(k,:) = idx';
end
end
