function [eig_bounds,U,K_CG,M_CG,A_xx,A_xy,A_yy, ...
    U_with_bdry,meshCG,diagnostics] = ...
    calc_eigen_bounds_cluster_robust( ...
    neig,triangle,N_LG,N_rho,LagrangeOrder)
%CALC_EIGEN_BOUNDS_CLUSTER_ROBUST  Adjacent bounds without an LG shift.
%
% For k=1,...,NEIG, the Crouzeix--Raviart/Liu lower bound
%
%   lambda_k >= lambda_CR,k/(1+C_h^2 lambda_CR,k)
%
% is paired with the conforming Lagrange Rayleigh upper bound.  Unlike a
% Lehmann--Goerisch shift above lambda_NEIG, this construction does not
% require separation from lambda_{NEIG+1}.  It therefore remains valid
% when the next eigenvalue is multiple or its cluster dimension is unknown.

global INTERVAL_MODE
if isempty(INTERVAL_MODE)
    INTERVAL_MODE = 0;
end

triangle = I_intval(triangle);
a = triangle(5);
b = triangle(6);

meshCR = make_mesh_by_gmsh(a,b,1/N_rho);
vertCR = meshCR.nodes;
edgeCR = meshCR.edges;
triCR = meshCR.elements;
bdCR = meshCR.boundary_edges;

isBnd = ismember(edgeCR,bdCR,'rows');
triByEdge = find_tri2edge(triCR,edgeCR);
bdEdgeIds = find(isBnd > 0);
[M_CR,K_CR] = create_matrix_crouzeix_raviart( ...
    triCR,edgeCR,vertCR,triByEdge);
dof = 1:size(edgeCR,1);
dof(bdEdgeIds) = [];
M_CR = M_CR(dof,dof);
K_CR = K_CR(dof,dof);

hmax = find_mesh_hmax(vertCR,edgeCR);
if INTERVAL_MODE
    [lambda_CR,lambda_CR_indices] = ...
        veigs(K_CR,M_CR,neig,'sm');
    lambda_CR = lambda_CR(:);
    lambda_CR_indices = lambda_CR_indices(:);
    % INTLAB overloads NUMEL(intval) as one even for a non-scalar
    % enclosure.  Count the ordinary midpoint array instead.
    n_cr = numel(I_mid(lambda_CR));
    if n_cr ~= numel(lambda_CR_indices) ...
            || ~all(ismember((1:neig)',lambda_CR_indices))
        error('calc_eigen_bounds_cluster_robust:BadCREigenvalueIndices', ...
            ['VEIGS did not return enclosures for precisely the ', ...
             'required low CR spectral indices.']);
    end
    ordered_CR = I_intval(zeros(neig,1));
    for k = 1:neig
        position = find(lambda_CR_indices == k,1);
        ordered_CR(k) = lambda_CR(position);
    end
    lambda_CR = ordered_CR;
    lambda_CR_indices = (1:neig)';
    Ch = I_intval('0.1893')*hmax;
else
    lambda_CR = eigs(K_CR,M_CR,neig,'sm');
    [~,order_CR] = sort(real(lambda_CR),'ascend');
    lambda_CR = lambda_CR(order_CR);
    lambda_CR_indices = (1:neig)';
    Ch = I_intval('0.1893')*hmax;
end
liu_candidates = lambda_CR./(1+lambda_CR.*Ch^2);
if numel(I_mid(liu_candidates)) < neig
    error('calc_eigen_bounds_cluster_robust:TooFewCREigenvalues', ...
        'The CR solve returned too few eigenvalues.');
end
liu_lower = zeros(neig,1);
for k = 1:neig
    liu_lower(k) = I_inf(liu_candidates(k));
end

meshCG = make_mesh_by_gmsh(a,b,1/N_LG);
[lambda_CG,U,U_with_bdry,K_CG,M_CG,A_xx,A_xy,A_yy,~] = ...
    laplace_eig_lagrange_detailed( ...
        LagrangeOrder,meshCG.nodes,meshCG.edges, ...
        meshCG.elements,meshCG.boundary_edges,neig);
lambda_CG = lambda_CG(:);
ritz_CG = verified_ritz_enclosures(U,K_CG,M_CG,neig);
lambda_CG_upper = zeros(neig,1);
for k = 1:neig
    lambda_CG_upper(k) = I_sup(ritz_CG(k));
end

eig_bounds = I_intval(zeros(neig,1));
for k = 1:neig
    lower = liu_lower(k);
    upper = lambda_CG_upper(k);
    if lower > upper
        error('calc_eigen_bounds_cluster_robust:InconsistentBounds', ...
            ['The Liu lower bound for lambda_%d exceeds its ', ...
             'conforming upper bound.'],k);
    end
    eig_bounds(k) = I_infsup(lower,upper);
end

diagnostics = struct();
diagnostics.method = 'CR-Liu-lower-plus-conforming-upper';
diagnostics.neig = neig;
diagnostics.lambda_CR = lambda_CR;
diagnostics.lambda_CR_indices = lambda_CR_indices;
diagnostics.Ch = Ch;
diagnostics.liu_lower_candidates = liu_candidates;
diagnostics.liu_lower_sorted = liu_lower;
diagnostics.lambda_CG = lambda_CG;
diagnostics.verified_ritz_CG = ritz_CG;
diagnostics.lambda_CG_upper_indexed = lambda_CG_upper;
% Backward-compatible alias retained for earlier diagnostics consumers.
diagnostics.lambda_CG_upper_sorted = lambda_CG_upper;
diagnostics.requires_upper_cluster_dimension = false;
diagnostics.requires_Lehmann_Goerisch_shift = false;
end
