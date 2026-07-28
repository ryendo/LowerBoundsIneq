function [eig_bounds,diagnostics] = cell_lower_eig_bound(cell_data)

    diagnostics = struct();
    cell_started = tic;

    % The eigenvalue lower bound is evaluated at the outer cell corner.
    x2 = I_intval(cell_data.x_sup);
    t2 = I_intval(cell_data.theta_sup);
    
    a4 = x2;  b4 = x2 * tan(t2);
    
    % Extract mesh/FEM parameters    
    mesh_size_rho = cell_data.mesh_size_lower_cr;
    
    isLG = I_mid(cell_data.isLG);

    if isLG
        % The Lehmann--Goerisch shift is always the validated CR--Liu
        % lower bound for lambda_2.
        [eig_bounds,diagnostics] = local_lg_lower_bound( ...
            a4,b4,mesh_size_rho,cell_data);

    else    
        % =========================================================
        % Standard CR Mode (isLG == 0)
        % =========================================================
        % Omega_mid consumes only lambda_1.  VEIGS is still asked for two
        % values, which supplies the adjacent spectral information used by
        % its verified local inclusion while avoiding the legacy four-value
        % solve.
        neig = 1;
        a_ = a4;
        b_ = b4;
        
        mesh_rho = make_mesh_by_gmsh(a_, b_, mesh_size_rho);
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
            [CR_eig,CR_indices] = ...
                veigs(CR_A1,CR_A0,neig+1,'sm');
            CR_eig = CR_eig(:);
            CR_indices = CR_indices(:);
            if numel(I_mid(CR_eig)) ~= numel(CR_indices) ...
                    || ~all(ismember((1:(neig+1))',CR_indices))
                error('cell_lower_eig_bound:BadCRIndices', ...
                    'The verified CR solve did not return indices 1:%d.', ...
                    neig+1);
            end
            ordered_CR = I_intval(zeros(neig+1,1));
            for k = 1:neig+1
                ordered_CR(k) = CR_eig(find(CR_indices==k,1));
            end
            CR_eig = ordered_CR;
        else
            CR_eig = eigs(CR_A1,CR_A0,neig+1,'sm');
            [~,order_CR] = sort(real(CR_eig),'ascend');
            CR_eig = CR_eig(order_CR);
        end
        
        Ch_val = I_intval('0.1893') * hmax;
        
        % Compute validated lower bounds using CR error estimate
        eig_bounds_full = CR_eig ./ (1 + CR_eig .* (Ch_val^2));

        eig_bounds = I_inf(eig_bounds_full(1));
        if ~isfinite(eig_bounds) || ~(eig_bounds > 0)
            error('cell_lower_eig_bound:InvalidCRBound', ...
                'The CR lower bound is not finite positive.');
        end
        diagnostics.method = 'CR';
        diagnostics.rho = I_intval(NaN);
        diagnostics.verified_ritz_upper = I_intval(NaN);
        diagnostics.separation_margin = I_intval(Inf);
        diagnostics.mu_upper = I_intval(-Inf);
        diagnostics.lower_bound = eig_bounds;
        diagnostics.shift_method = 'not-applicable';
        diagnostics.fallback_used = false;
        diagnostics.fallback_reason = '';
        diagnostics.cr_shift = I_intval(NaN);
        diagnostics.B_verified_spd = true;
        diagnostics.timing = struct( ...
            'lg_setup_seconds',0, ...
            'cr_shift_seconds',0, ...
            'lg_solve_seconds',0, ...
            'total_seconds',toc(cell_started));
        diagnostics.conditions_verified = true;
    end
end


function [eig_bound,diagnostics] = local_lg_lower_bound( ...
    a,b,mesh_size_rho,cell_data)
total_started = tic;

if ~isfield(cell_data,'mesh_size_lower_LG') ...
        || ~isfield(cell_data,'fem_order_lower_LG')
    error('cell_lower_eig_bound:MissingLGParameters', ...
        'The LG branch requires mesh_size_lower_LG and fem_order_lower_LG.');
end
mesh_size_LG = cell_data.mesh_size_lower_LG;
Lagrange_order = cell_data.fem_order_lower_LG;

cr_started = tic;
[rho,cr_shift_info] = ...
    local_cr_lambda2_shift(a,b,mesh_size_rho);
cr_shift_seconds = toc(cr_started);

setup_started = tic;
mesh_LG = make_mesh_by_gmsh(a,b,mesh_size_LG);
[~,LA_eigf,LA_eigf_with_bdry,LA_A,LA_M,~,~,~,~] = ...
    laplace_eig_lagrange_detailed( ...
        Lagrange_order,mesh_LG.nodes,mesh_LG.edges, ...
        mesh_LG.elements,mesh_LG.boundary_edges,1);
ritz = verified_ritz_enclosures(LA_eigf,LA_A,LA_M,1);
Lambda = I_intval(I_sup(ritz(1)));
A0 = LA_eigf'*LA_A*LA_eigf;
A1 = LA_eigf'*LA_M*LA_eigf;
A2 = RT_Hdiv_problem_dirichlet( ...
    mesh_LG,Lagrange_order,LA_eigf_with_bdry);
lg_setup_seconds = toc(setup_started);

solve_started = tic;
[eig_bound,final_certificate] = ...
    local_apply_lg_shift(A0,A1,A2,rho,Lambda);
lg_solve_seconds = toc(solve_started);

diagnostics = struct();
diagnostics.method = 'LG';
diagnostics.shift_method = 'CR-Liu';
diagnostics.rho = rho;
diagnostics.verified_ritz_upper = Lambda;
diagnostics.separation_margin = ...
    final_certificate.separation_margin;
diagnostics.mu_upper = final_certificate.mu_upper;
diagnostics.B_verified_spd = ...
    logical(final_certificate.B_verified_spd);
diagnostics.lower_bound = eig_bound;
diagnostics.conditions_verified = true;
diagnostics.fallback_used = false;
diagnostics.fallback_reason = '';
diagnostics.cr_shift = rho;
diagnostics.cr_shift_info = cr_shift_info;
diagnostics.timing = struct( ...
    'lg_setup_seconds',lg_setup_seconds, ...
    'cr_shift_seconds',cr_shift_seconds, ...
    'lg_solve_seconds',lg_solve_seconds, ...
    'total_seconds',toc(total_started));
end


function [eig_bound,certificate] = ...
    local_apply_lg_shift(A0,A1,A2,rho,Lambda)
global INTERVAL_MODE
rho = I_intval(I_inf(rho));
Lambda = I_intval(I_sup(Lambda));
if ~(I_sup(Lambda) < I_inf(rho))
    error('cell_lower_eig_bound:SeparationFailed', ...
        ['LG separation failed: the verified conforming Ritz ', ...
         'upper bound must be below the lambda_2 lower shift.']);
end

AL = I_hull(A0-rho*A1,(A0-rho*A1)');
BL_raw = A0-2*rho*A1+rho*rho*A2;
BL = I_hull(BL_raw,BL_raw');
if INTERVAL_MODE && ~isspd(BL)
    error('cell_lower_eig_bound:BLNotSPD', ...
        'The LG matrix B was not verified positive definite.');
end

mus = I_eig(AL,BL,1);
if any(I_sup(mus) >= 0)
    error('cell_lower_eig_bound:MuNotNegative', ...
        'The LG generalized eigenvalue was not verified negative.');
end
lg_lower = rho-rho./(1-mus(end:-1:1));
eig_bound = I_inf(lg_lower(1));
if ~isfinite(eig_bound) || ~(eig_bound > 0)
    error('cell_lower_eig_bound:InvalidLGBound', ...
        'The LG lower bound is not finite positive.');
end

certificate = struct();
certificate.separation_margin = ...
    I_intval(I_inf(rho-Lambda));
certificate.mu_upper = I_intval(max(I_sup(mus)));
certificate.B_verified_spd = ...
    logical(~INTERVAL_MODE || isspd(BL));
end


function [rho,info] = local_cr_lambda2_shift(a,b,mesh_size)
global INTERVAL_MODE
mesh = make_mesh_by_gmsh(a,b,mesh_size);
is_boundary = ismember( ...
    mesh.edges,mesh.boundary_edges,'rows');
tri_by_edge = find_tri2edge( ...
    mesh.elements,mesh.edges);
[Mcr,Kcr] = create_matrix_crouzeix_raviart( ...
    mesh.elements,mesh.edges,mesh.nodes,tri_by_edge);
dofs = 1:size(mesh.edges,1);
dofs(is_boundary > 0) = [];
Mcr = Mcr(dofs,dofs);
Kcr = Kcr(dofs,dofs);
hmax = find_mesh_hmax(mesh.nodes,mesh.edges);

if INTERVAL_MODE
    [values,indices] = veigs(Kcr,Mcr,2,'sm');
    values = values(:);
    indices = indices(:);
    if numel(I_mid(values)) ~= numel(indices)
        error('cell_lower_eig_bound:BadCRIndices', ...
            'The verified CR solve returned incoherent indices.');
    end
    ordered = I_intval(zeros(2,1));
    for k = 1:2
        location = find(indices == k);
        if numel(location) ~= 1
            error('cell_lower_eig_bound:BadCRIndices', ...
                'The verified CR solve did not return spectral index %d.',k);
        end
        ordered(k) = values(location);
    end
    values = ordered;
else
    values = eigs(Kcr,Mcr,2,'sm');
    [~,order] = sort(real(values),'ascend');
    values = values(order);
end

Ch = I_intval('0.1893')*hmax;
lower_bounds = values./(1+values*(Ch^2));
rho = I_intval(I_inf(lower_bounds(2)));
if ~(I_inf(rho) > 0) || ~isfinite(I_sup(rho))
    error('cell_lower_eig_bound:BadShift', ...
        'The CR lower shift for lambda_2 is not finite positive.');
end
info = struct( ...
    'used',true, ...
    'method','Crouzeix-Raviart-Liu', ...
    'hmax',hmax, ...
    'lambda2_discrete_enclosure',values(2), ...
    'rho',rho);
end
