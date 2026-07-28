function [residual_norm_upper, diagnostics] = ...
    RT_shifted_equilibrated_flux_dirichlet(mesh, RT_order, ...
    Lagrange_order, u_full, v_full, P_e, alpha_h, lambda_h,options)
%RT_SHIFTED_EQUILIBRATED_FLUX_DIRICHLET  Shifted RT flux majorant.
%
% Find sigma_h in RT_k which minimizes
%
%   ||sigma_h-(P_e grad u_h+grad v_h)||_{L2}
%
% subject to
%
%   div sigma_h = -(alpha_h u_h+lambda_h v_h).
%
% The scalar source is represented exactly by div(RT_k), so this routine
% currently requires k to equal the Lagrange degree.
%
% The default midpoint-defect strategy fixes the midpoint mixed solution
% sigma_0 and rigorously bounds both
%
%   eta_grad = ||sigma_0-(P_e grad u_h+grad v_h)||,
%   eta_src  = ||div sigma_0+alpha_h u_h+lambda_h v_h||.
%
% Given OPTIONS.lambda1_lower, the first output is the energy-dual bound
%
%   eta_grad + eta_src/sqrt(lambda1_lower).
%
% OPTIONS.solve_strategy='verified-schur' retains the exact-equilibrium
% interval solve used by the original implementation.
%
% u_full and v_full contain all Lagrange degrees of freedom, including
% zeros on the Dirichlet boundary.  P_e is the symmetric matrix in the
% first shape derivative form a_e(w,z)=(P_e grad w,grad z).

global INTERVAL_MODE
if isempty(INTERVAL_MODE)
    INTERVAL_MODE = 0;
end
if nargin < 9 || isempty(options)
    options = struct();
end
if isfield(options,'solve_strategy') && ~isempty(options.solve_strategy)
    requested_strategy = lower(char(options.solve_strategy));
else
    requested_strategy = 'midpoint-defect';
end
if any(strcmp(requested_strategy,{'verified-schur-complement','exact'}))
    requested_strategy = 'verified-schur';
end
if ~any(strcmp(requested_strategy,{'midpoint-defect','verified-schur'}))
    error('RT_shifted_equilibrated_flux_dirichlet:BadSolveStrategy', ...
        ['options.solve_strategy must be ''midpoint-defect'' or ', ...
         '''verified-schur''.']);
end
has_lambda1_lower = isfield(options,'lambda1_lower') ...
    && ~isempty(options.lambda1_lower);

if RT_order ~= Lagrange_order
    error('RT_shifted_equilibrated_flux_dirichlet:OrderMismatch', ...
        ['Exact equilibration of a degree-p Lagrange source currently ', ...
         'requires RT_order=Lagrange_order.']);
end
if RT_order < 1
    error('RT_shifted_equilibrated_flux_dirichlet:OrderTooLow', ...
        'The material-derivative source requires order at least one.');
end

vert = I_intval(mesh.nodes);
edge = mesh.edges;
tri = mesh.elements;
nt = size(tri, 1);
ne = size(edge, 1);
[tri, tri2edge] = local_find_tri2edge(edge, tri, ne, nt);

ndof_lg = size(u_full, 1);
if size(v_full, 1) ~= ndof_lg
    error('RT_shifted_equilibrated_flux_dirichlet:SizeMismatch', ...
        'u_full and v_full must have the same number of entries.');
end

[A, Bdiv, target_load, source_load, target_norm_sq, ...
    Wdg, assembly] = ...
    assemble_shifted_problem(vert, edge, tri, tri2edge, RT_order, ...
    u_full(:), v_full(:), P_e, alpha_h, lambda_h);

n_sigma = size(A, 1);
n_dg = size(Bdiv, 2);
if strcmp(requested_strategy,'verified-schur') && INTERVAL_MODE
    % The Schur form follows the verified path already used by
    % RT_Hdiv_problem_dirichlet and avoids a large indefinite interval
    % inverse.  From
    %   A*sigma+B*y=H,  B'*sigma=-F
    % we obtain
    %   (B'*A^{-1}B)y=F+B'*A^{-1}H.
    AinvB = I_solve(A,Bdiv);
    AinvH = I_solve(A,target_load);
    schur = Bdiv'*AinvB;
    multiplier = I_solve(schur,source_load+Bdiv'*AinvH);
    sigma = AinvH-AinvB*multiplier;
    solve_strategy = 'verified-Schur-complement';
elseif strcmp(requested_strategy,'verified-schur')
    mixed = [A, Bdiv; Bdiv', I_zeros(n_dg, n_dg)];
    rhs = [target_load; -source_load];
    solution = I_solve(mixed,rhs);
    sigma = solution(1:n_sigma);
    multiplier = solution(n_sigma+1:end);
    solve_strategy = 'direct-mixed-double';
else
    % A fixed double vector has no interval solve-box wrapping.  Its
    % nonzero equilibrium defect is measured rigorously below.
    A0 = I_mid(A);
    B0 = I_mid(Bdiv);
    H0 = I_mid(target_load);
    F0 = I_mid(source_load);
    mixed0 = [A0, B0; B0', sparse(n_dg,n_dg)];
    rhs0 = [H0; -F0];
    solution0 = mixed0\rhs0;
    sigma = I_intval(solution0(1:n_sigma));
    multiplier = I_intval(solution0(n_sigma+1:end));
    solve_strategy = 'midpoint-mixed-with-verified-defect';
end

% In the verified-Schur branch the box contains an exactly equilibrated
% solution.  In the midpoint-defect branch SIGMA is simply a fixed trial
% vector; its quadratic gradient mismatch is still evaluated rigorously,
% and its nonzero equilibrium defect is added below.
residual_sq_enclosure = sigma' * A * sigma ...
    - 2 * sigma' * target_load + target_norm_sq;
residual_sq_upper = max(I_sup(residual_sq_enclosure), 0);
if isnan(residual_sq_upper) || isinf(residual_sq_upper)
    error('RT_shifted_equilibrated_flux_dirichlet:NonfiniteResidual', ...
        'The verified mixed solve did not produce a finite flux majorant.');
end
gradient_residual = I_intval( ...
    I_sup(sqrt(I_intval(residual_sq_upper))));

equilibrium_defect = Bdiv' * sigma + source_load;
stationarity_defect = A * sigma + Bdiv * multiplier - target_load;

if strcmp(requested_strategy,'midpoint-defect')
    source_residual = local_DG_dual_mass_norm_upper( ...
        Wdg,equilibrium_defect,assembly.local_Lagrange_basis_size,nt);
else
    % The exact solution enclosed by the verified Schur box satisfies the
    % equilibrium equation exactly; the displayed defect box is only an
    % enclosure-evaluation diagnostic.
    source_residual = I_intval(0);
end

if has_lambda1_lower
    lambda1_lower_input = I_intval(options.lambda1_lower);
    lambda1_lower = I_intval(I_inf(lambda1_lower_input));
    if ~(I_inf(lambda1_lower) > 0)
        error('RT_shifted_equilibrated_flux_dirichlet:BadPoincareConstant', ...
            'options.lambda1_lower must be strictly positive.');
    end
    source_dual_residual = source_residual/sqrt(lambda1_lower);
    dual_residual = I_intval(I_sup( ...
        gradient_residual+source_dual_residual));
elseif I_sup(source_residual) == 0
    lambda1_lower = I_intval(NaN);
    source_dual_residual = I_intval(0);
    dual_residual = gradient_residual;
else
    error('RT_shifted_equilibrated_flux_dirichlet:MissingPoincareConstant', ...
        ['The midpoint-defect strategy requires ', ...
         'options.lambda1_lower.']);
end
residual_norm_upper = dual_residual;

diagnostics = struct();
diagnostics.rigorous = logical(INTERVAL_MODE);
diagnostics.RT_order = RT_order;
diagnostics.Lagrange_order = Lagrange_order;
diagnostics.num_RT_dofs = n_sigma;
diagnostics.num_DG_dofs = n_dg;
diagnostics.num_elements = nt;
diagnostics.solve_strategy = solve_strategy;
diagnostics.requested_solve_strategy = requested_strategy;
diagnostics.target_norm_sq = target_norm_sq;
diagnostics.residual_sq_enclosure = residual_sq_enclosure;
diagnostics.residual_sq_upper = I_intval(residual_sq_upper);
diagnostics.gradient_residual = gradient_residual;
diagnostics.source_residual = source_residual;
diagnostics.lambda1_lower = lambda1_lower;
diagnostics.source_dual_residual = source_dual_residual;
diagnostics.dual_residual = dual_residual;
diagnostics.residual_norm_upper = residual_norm_upper;
diagnostics.equilibrium_defect_norm_upper = ...
    interval_vector_norm_upper(equilibrium_defect);
diagnostics.stationarity_defect_norm_upper = ...
    interval_vector_norm_upper(stationarity_defect);
diagnostics.sigma = sigma;
diagnostics.multiplier = multiplier;
diagnostics.assembly = assembly;
end


function [Aglob, Bglob, Hglob, Fglob, Gnorm2, Wglob, diagnostics] = ...
    assemble_shifted_problem(vert, edge, tri, tri2edge, order, ...
    u_full, v_full, P_e, alpha_h, lambda_h)

[rt_abc, rt_ijk, nrt] = local_RT_basis(order);
[lg_basis, nlg] = local_Lagrange_basis(order);
nmono_rt = local_poly_nbasis(order + 1);
nmono_lg = local_poly_nbasis(order);

M_rt = local_cross_monomial_mass(order + 1, order + 1);
M_lg = local_cross_monomial_mass(order, order);
M_cross = local_cross_monomial_mass(order + 1, order);

rt_value = cell(nrt, 1);
rt_div = cell(nrt, 1);
for a = 1:nrt
    rt_value{a} = local_RT_coord_basis(rt_abc, rt_ijk, a, order);
    rt_div{a} = local_RT_div_basis(rt_abc, rt_ijk, a, order);
end

lg_value = cell(nlg, 1);
lg_grad = cell(nlg, 1);
for a = 1:nlg
    lg_value{a} = local_Lagrange_coord_basis(lg_basis, a, order);
    lg_grad{a} = local_Lagrange_grad_basis(lg_basis, a, order);
end

Aref_tensor = cell(nrt, nrt);
Bref = I_zeros(nrt, nlg);
Mlg_basis = I_zeros(nlg, nlg);
for a = 1:nrt
    for b = 1:a
        Aref_tensor{a,b} = rt_value{a}' * M_rt * rt_value{b};
    end
    for b = 1:nlg
        Bref(a,b) = rt_div{a}' * M_lg * lg_value{b};
    end
end
for a = 1:nlg
    for b = 1:nlg
        Mlg_basis(a,b) = lg_value{a}' * M_lg * lg_value{b};
    end
end

ne = size(edge, 1);
nv = size(vert, 1);
nt = size(tri, 1);
n_rt_global = (order+1)*ne + order*(order+1)*nt;
n_dg_global = nlg*nt;

Aglob = I_sparse(n_rt_global, n_rt_global);
Bglob = I_sparse(n_rt_global, n_dg_global);
Hglob = I_zeros(n_rt_global, 1);
Fglob = I_zeros(n_dg_global, 1);
Wglob = I_sparse(n_dg_global,n_dg_global);
Gnorm2 = I_intval(0);

min_det_lower = Inf;
max_det_upper = 0;

for elem = 1:nt
    verts = tri(elem,:);
    x1 = vert(verts(1),1); y1 = vert(verts(1),2);
    x2 = vert(verts(2),1); y2 = vert(verts(2),2);
    x3 = vert(verts(3),1); y3 = vert(verts(3),2);
    Bmap = [x2-x1, x3-x1; y2-y1, y3-y1];
    detB = det(Bmap);
    if ~(I_inf(detB) > 0)
        error('RT_shifted_equilibrated_flux_dirichlet:BadElement', ...
            'Element %d is not certified counter-clockwise and nondegenerate.', elem);
    end
    min_det_lower = min(min_det_lower, I_inf(detB));
    max_det_upper = max(max_det_upper, I_sup(detB));
    Binv = [y3-y1, x1-x3; y1-y2, x2-x1] / detB;

    [rt_map, rt_sign] = local_RT_dof_map( ...
        order, edge, tri2edge(elem,:), verts, elem, ne);
    lg_map = local_Lagrange_dof_map( ...
        order, edge, tri2edge(elem,:), verts, elem, nv, ne);
    dg_map = (elem-1)*nlg + (1:nlg);

    if max(lg_map) > size(u_full,1)
        error('RT_shifted_equilibrated_flux_dirichlet:BadLagrangeVector', ...
            'Full Lagrange vector is shorter than the mesh DOF map.');
    end

    ucoef = u_full(lg_map);
    vcoef = v_full(lg_map);

    grad_u_ref = I_zeros(nmono_lg, 2);
    grad_v_ref = I_zeros(nmono_lg, 2);
    for a = 1:nlg
        grad_u_ref = grad_u_ref + lg_grad{a} * ucoef(a);
        grad_v_ref = grad_v_ref + lg_grad{a} * vcoef(a);
    end

    % Row-vector convention: grad_phys^T=grad_ref^T*B^{-1}.
    Gphys = grad_u_ref * Binv * P_e' + grad_v_ref * Binv;
    % The Piola pairing cancels det(B):
    % int (B tauhat/detB).Gphys detB = int tauhat.(B'Gphys).
    Gpiola_pair = Gphys * Bmap;

    Alocal = I_zeros(nrt, nrt);
    Hlocal = I_zeros(nrt, 1);
    for a = 1:nrt
        for b = 1:a
            Alocal(a,b) = trace(Bmap * Aref_tensor{a,b} * Bmap') / detB;
        end
        Hlocal(a) = trace(rt_value{a}' * M_cross * Gpiola_pair);
    end
    Alocal = Alocal + tril(Alocal,-1)';

    fcoef = alpha_h * ucoef + lambda_h * vcoef;
    Wlocal = Mlg_basis*detB;
    Flocal = Wlocal*fcoef;
    Gnorm2 = Gnorm2 + trace(Gphys' * M_lg * Gphys) * detB;

    S = diag(rt_sign);
    Aglob(rt_map,rt_map) = Aglob(rt_map,rt_map) + S*Alocal*S;
    Bglob(rt_map,dg_map) = Bglob(rt_map,dg_map) + S*Bref;
    Hglob(rt_map) = Hglob(rt_map) + S*Hlocal;
    Fglob(dg_map) = Fglob(dg_map) + Flocal;
    Wglob(dg_map,dg_map) = Wlocal;
end

diagnostics = struct();
diagnostics.local_RT_basis_size = nrt;
diagnostics.local_Lagrange_basis_size = nlg;
diagnostics.nmono_RT = nmono_rt;
diagnostics.nmono_Lagrange = nmono_lg;
diagnostics.min_element_det_lower = min_det_lower;
diagnostics.max_element_det_upper = max_det_upper;
end


function [map, signs] = local_RT_dof_map(order, edge, elem_edges, ...
    verts, elem, ne)
nrt = local_RT_nbasis(order);
map = zeros(nrt,1);
signs = ones(nrt,1);
local_edge_start_vertex = [2,3,1];
for side = 1:3
    rows = (order+1)*(side-1)+(1:(order+1));
    edge_id = elem_edges(side);
    if edge(edge_id,1) == verts(local_edge_start_vertex(side))
        map(rows) = (order+1)*(edge_id-1)+(1:(order+1));
    else
        map(rows) = (order+1)*edge_id:-1:(order+1)*(edge_id-1)+1;
        signs(rows) = -1;
    end
end
interior = (order+1)*3+1:nrt;
map(interior) = (order+1)*ne ...
    + order*(order+1)*(elem-1)+(1:(order*(order+1)));
end


function map = local_Lagrange_dof_map(order, edge, elem_edges, ...
    verts, elem, nv, ne)
nlg = local_poly_nbasis(order);
map = zeros(nlg,1);
map(1:3) = verts;
local_edge_start_vertex = [2,3,1];
for side = 1:3
    rows = 3+(order-1)*(side-1)+(1:(order-1));
    if isempty(rows)
        continue;
    end
    edge_id = elem_edges(side);
    first = nv+(order-1)*(edge_id-1)+1;
    last = nv+(order-1)*edge_id;
    if edge(edge_id,1) == verts(local_edge_start_vertex(side))
        map(rows) = first:last;
    else
        map(rows) = last:-1:first;
    end
end
ninterior = (order-1)*(order-2)/2;
if ninterior > 0
    rows = 3+3*(order-1)+(1:ninterior);
    map(rows) = nv+(order-1)*ne+ninterior*(elem-1)+(1:ninterior);
end
end


function [tri, tri2edge] = local_find_tri2edge(edge, tri, ne, nt)
tri2edge = zeros(nt,3);
edge_key = sort(edge,2)*[ne;1];
local_edges = [2 3;1 3;1 2];
for elem = 1:nt
    value = sort(reshape(tri(elem,local_edges),3,2),2)*[ne;1];
    [found, idx] = ismember(value,edge_key);
    if ~all(found)
        error('RT_shifted_equilibrated_flux_dirichlet:MissingEdge', ...
            'The edge list does not contain all element edges.');
    end
    tri2edge(elem,:) = idx';
end
end


function M = local_cross_monomial_mass(degree_left, degree_right)
left = local_create_ijk(degree_left);
right = local_create_ijk(degree_right);
M = I_zeros(size(left,1),size(right,1));
for a = 1:size(left,1)
    for b = 1:size(right,1)
        exponent = left(a,:)+right(b,:);
        numerator = factorial(exponent(1))*factorial(exponent(2)) ...
            * factorial(exponent(3));
        denominator = factorial(sum(exponent)+2);
        M(a,b) = I_intval(numerator)/I_intval(denominator);
    end
end
end


function e = local_RT_coord_basis(abc, ijk, idx, order)
e = I_zeros(local_poly_nbasis(order+1),2);
a = abc(idx,1); b = abc(idx,2); c = abc(idx,3);
i = ijk(idx,1); j = ijk(idx,2); k = ijk(idx,3);
e(local_map_ijk(i+1,j,k,order+1),1) = a;
e(local_map_ijk(i,j+1,k,order+1),1) = a+c;
e(local_map_ijk(i,j,k+1,order+1),1) = a;
e(local_map_ijk(i+1,j,k,order+1),2) = b;
e(local_map_ijk(i,j+1,k,order+1),2) = b;
e(local_map_ijk(i,j,k+1,order+1),2) = b+c;
end


function e = local_RT_div_basis(abc, ijk, idx, order)
e = I_zeros(local_poly_nbasis(order),1);
f = I_zeros(local_poly_nbasis(order),1);
a = abc(idx,1); b = abc(idx,2); c = abc(idx,3);
i = ijk(idx,1); j = ijk(idx,2); k = ijk(idx,3);
e(local_map_ijk(i,j,k,order)) = a*(j-i)+c*(j+1);
e(local_map_ijk(i-1,j+1,k,order)) = -(a+c)*i;
e(local_map_ijk(i-1,j,k+1,order)) = -a*i;
e(local_map_ijk(i+1,j-1,k,order)) = a*j;
e(local_map_ijk(i,j-1,k+1,order)) = a*j;
f(local_map_ijk(i,j,k,order)) = b*(k-i)+c*(k+1);
f(local_map_ijk(i-1,j,k+1,order)) = -(b+c)*i;
f(local_map_ijk(i-1,j+1,k,order)) = -b*i;
f(local_map_ijk(i+1,j,k-1,order)) = b*k;
f(local_map_ijk(i,j+1,k-1,order)) = b*k;
e = e+f;
end


function e = local_Lagrange_coord_basis(basis, idx, order)
e = I_zeros(local_poly_nbasis(order),1);
i = basis(idx,1); j = basis(idx,2); k = basis(idx,3);
e(local_map_ijk(i,j,k,order)) = 1;
end


function e = local_Lagrange_grad_basis(basis, idx, order)
n = local_poly_nbasis(order);
dx = I_zeros(n,1);
dy = I_zeros(n,1);
i = basis(idx,1); j = basis(idx,2); k = basis(idx,3);
dx(local_map_ijk(i,j,k,order)) = -i+j;
dx(local_map_ijk(i-1,j+1,k,order)) = -i;
dx(local_map_ijk(i-1,j,k+1,order)) = -i;
dx(local_map_ijk(i+1,j-1,k,order)) = j;
dx(local_map_ijk(i,j-1,k+1,order)) = j;
dy(local_map_ijk(i,j,k,order)) = -i+k;
dy(local_map_ijk(i-1,j+1,k,order)) = -i;
dy(local_map_ijk(i-1,j,k+1,order)) = -i;
dy(local_map_ijk(i+1,j,k-1,order)) = k;
dy(local_map_ijk(i,j+1,k-1,order)) = k;
e = [dx,dy];
end


function idx = local_map_ijk(i,j,k,degree)
if i < 0 || j < 0 || k < 0
    idx = [];
    return;
end
if i+j+k ~= degree
    idx = [];
    return;
end
idx = (degree-i)*(degree-i+1)/2+(degree-i-j)+1;
end


function ijk = local_create_ijk(degree)
ijk = zeros(local_poly_nbasis(degree),3);
idx = 1;
for i = degree:-1:0
    for j = degree-i:-1:0
        ijk(idx,:) = [i,j,degree-i-j];
        idx = idx+1;
    end
end
end


function [basis, nbasis] = local_Lagrange_basis(order)
nbasis = local_poly_nbasis(order);
basis = zeros(nbasis,3);
idx = 1;
basis(idx,:) = [order,0,0]; idx = idx+1;
if order > 0
    basis(idx,:) = [0,order,0]; idx = idx+1;
    basis(idx,:) = [0,0,order]; idx = idx+1;
end
for p = order-1:-1:1
    basis(idx,:) = [0,p,order-p]; idx = idx+1;
end
for p = order-1:-1:1
    basis(idx,:) = [order-p,0,p]; idx = idx+1;
end
for p = order-1:-1:1
    basis(idx,:) = [p,order-p,0]; idx = idx+1;
end
for p = order-2:-1:1
    for q = order-1-p:-1:1
        basis(idx,:) = [p,q,order-p-q]; idx = idx+1;
    end
end
end


function [abc, ijk, nbasis] = local_RT_basis(order)
nbasis = local_RT_nbasis(order);
abc = zeros(nbasis,3);
ijk = zeros(nbasis,3);
idx = 1;
abc(idx,:)=[0,0,1]; ijk(idx,:)=[0,order,0]; idx=idx+1;
for p=order-1:-1:1
    abc(idx,:)=[1,0,0]; ijk(idx,:)=[0,p,order-p]; idx=idx+1;
end
if order~=0
    abc(idx,:)=[0,1,0]; ijk(idx,:)=[0,0,order]; idx=idx+1;
end
abc(idx,:)=[-1,0,1]; ijk(idx,:)=[0,0,order]; idx=idx+1;
for p=order-1:-1:1
    abc(idx,:)=[-1,0,0]; ijk(idx,:)=[order-p,0,p]; idx=idx+1;
end
if order~=0
    abc(idx,:)=[-1,0,0]; ijk(idx,:)=[order,0,0]; idx=idx+1;
end
abc(idx,:)=[0,-1,1]; ijk(idx,:)=[order,0,0]; idx=idx+1;
for p=order-1:-1:1
    abc(idx,:)=[0,-1,0]; ijk(idx,:)=[p,order-p,0]; idx=idx+1;
end
if order~=0
    abc(idx,:)=[0,-1,1]; ijk(idx,:)=[0,order,0]; idx=idx+1;
end
if order>=1
    abc(idx,:)=[0,0,1]; ijk(idx,:)=[order,0,0]; idx=idx+1;
    abc(idx,:)=[-1,0,1]; ijk(idx,:)=[0,order,0]; idx=idx+1;
end
if order>=2
    for p=order-1:-1:1
        abc(idx,:)=[1,-1,0]; ijk(idx,:)=[0,p,order-p]; idx=idx+1;
    end
    for p=order-1:-1:1
        abc(idx,:)=[0,1,0]; ijk(idx,:)=[order-p,0,p]; idx=idx+1;
    end
    for p=order-1:-1:1
        abc(idx,:)=[-1,0,0]; ijk(idx,:)=[p,order-p,0]; idx=idx+1;
    end
    for p=order-1:-1:1
        abc(idx,:)=[0,0,1]; ijk(idx,:)=[p,order-p,0]; idx=idx+1;
    end
end
if order>=3
    for p=order-2:-1:1
        for q=order-1-p:-1:1
            abc(idx,:)=[1,0,0]; ijk(idx,:)=[p,q,order-p-q]; idx=idx+1;
        end
    end
    for p=order-2:-1:1
        for q=order-1-p:-1:1
            abc(idx,:)=[0,1,0]; ijk(idx,:)=[p,q,order-p-q]; idx=idx+1;
        end
    end
end
if idx-1 ~= nbasis
    error('RT_shifted_equilibrated_flux_dirichlet:InternalBasisError', ...
        'RT basis construction returned the wrong number of functions.');
end
end


function n = local_poly_nbasis(degree)
n = (degree+1)*(degree+2)/2;
end


function n = local_RT_nbasis(order)
n = (order+1)*(order+3);
end


function value = local_DG_dual_mass_norm_upper( ...
    Wglob,load_defect,nlocal,nelements)
% If d_j=(q,phi_j), then ||q||^2=d'W^{-1}d.  The discontinuous mass matrix
% is element block diagonal, so verified small-block solves avoid wrapping
% and the cost of a global interval inverse.
norm_sq = I_intval(0);
for elem = 1:nelements
    idx = (elem-1)*nlocal+(1:nlocal);
    Wlocal = Wglob(idx,idx);
    dlocal = load_defect(idx);
    qlocal = I_solve(Wlocal,dlocal);
    norm_sq = norm_sq+dlocal'*qlocal;
end
upper = max(I_sup(norm_sq),0);
if ~isfinite(upper)
    error('RT_shifted_equilibrated_flux_dirichlet:NonfiniteSourceDefect', ...
        'The verified DG source-defect norm is not finite.');
end
value = I_intval(I_sup(sqrt(I_intval(upper))));
end


function value = interval_vector_norm_upper(v)
if isempty(v)
    value = I_intval(0);
else
    value = I_intval(I_sup(norm(v,2)));
end
end
