function test_residual_hessian_estimator_double(requested_interval_mode)
%TEST_RESIDUAL_HESSIAN_ESTIMATOR_DOUBLE  Non-rigorous assembly sanity test.
%
% By default this test runs with INTERVAL_MODE=0.  Passing true also gives a
% small INTLAB assembly smoke test (the mock spectral data are still not a
% proof certificate).  It checks the bordered
% material-derivative solve, shifted RT equilibration, the ordering of both
% Hessian bounds, arbitrary-index gap constants (i=2), and the special
% C_minus=0 upper correction for i=1.  Certification is exercised on the
% HPC/INTLAB path; this local test is only meant to catch assembly and sign
% regressions without requiring INTLAB or gmsh.

this_file = mfilename('fullpath');
repo_root = fileparts(fileparts(this_file));
addpath(fullfile(repo_root,'src','interval'),'-begin');
addpath(fullfile(repo_root,'src','fem'),'-begin');
addpath(fullfile(repo_root,'src','lib','VFEM2D_revised'),'-begin');

global INTERVAL_MODE
if nargin < 1
    requested_interval_mode = false;
end
INTERVAL_MODE = logical(requested_interval_mode);

apex = [0.37,0.82];
mesh = make_uniform_triangle_mesh(apex,6);
p = 1;
neig = 3;
[lambda_h_all,U,U_full,K,M,Axx,Axy,Ayy] = ...
    laplace_eig_lagrange_detailed(p,mesh.nodes,mesh.edges, ...
        mesh.elements,mesh.boundary_edges,neig);

Pe = [0,-1/apex(2);-1/apex(2),0];
Pee = [2/apex(2)^2,0;0,0];
Ae = shape_matrix(Pe,Axx,Axy,Ayy);
Aee = shape_matrix(Pee,Axx,Axy,Ayy);

zero_full = I_zeros(size(U_full,1),1);
zero_flux_options = struct( ...
    'solve_strategy','midpoint-defect','lambda1_lower',1);
[zero_flux,~] = RT_shifted_equilibrated_flux_dirichlet( ...
    mesh,p,p,zero_full,zero_full,Pe,I_intval(0),I_intval(0), ...
    zero_flux_options);
assert(I_sup(zero_flux) == 0);

% Keep the original exact-equilibrium verified Schur path available as a
% regression/comparison option.
exact_flux_options = struct('solve_strategy','verified-schur');
[zero_flux_exact,exact_flux_info] = ...
    RT_shifted_equilibrated_flux_dirichlet( ...
        mesh,p,p,zero_full,zero_full,Pe,I_intval(0),I_intval(0), ...
        exact_flux_options);
assert(I_sup(zero_flux_exact) == 0);
assert(strcmp(exact_flux_info.requested_solve_strategy,'verified-schur'));

for i = 1:2
    u = U(:,i);
    scale = 1/sqrt(u'*M*u);
    u = u*scale;
    u_full = U_full(:,i)*scale;
    lambda_h = (u'*K*u)/(u'*M*u);
    alpha_h = (u'*Ae*u)/(u'*M*u);

    K0 = I_mid(K);
    M0 = I_mid(M);
    Ae0 = I_mid(Ae);
    u0 = I_mid(u);
    u0 = u0/sqrt(u0'*M0*u0);
    lambda0 = (u0'*K0*u0)/(u0'*M0*u0);
    alpha0 = (u0'*Ae0*u0)/(u0'*M0*u0);
    m_u = M0*u0;
    bordered = [K0-lambda0*M0,m_u;m_u',0];
    rhs = [-Ae0*u0+alpha0*m_u;0];
    sol = bordered\rhs;
    v = sol(1:end-1);
    assert(norm(bordered*sol-rhs,2) < 1e-8);

    boundary_vertices = unique(mesh.boundary_edges(:));
    inside = (1:size(mesh.nodes,1))';
    inside(boundary_vertices) = [];
    v_full = zeros(size(mesh.nodes,1),1);
    v_full(inside) = v;

    flux_options = struct( ...
        'solve_strategy','midpoint-defect', ...
        'lambda1_lower',I_inf(I_intval(0.9*lambda_h_all(1))));
    [eta_flux,flux_info] = ...
        RT_shifted_equilibrated_flux_dirichlet( ...
            mesh,p,p,u_full,v_full,Pe,alpha_h,lambda_h, ...
            flux_options);
    assert(isfinite(I_sup(eta_flux)) && I_inf(eta_flux) >= 0);
    assert(isfinite(I_sup(flux_info.residual_sq_upper)));

    if i == 1
        % Reuse of the same RT solver for the reference eigenproblem:
        % P_e=I, v=0, alpha=mu, lambda=0 gives
        % div(sigma)=-mu*u_h and targets grad(u_h).
        [rho_eig,eig_flux_info] = ...
            RT_shifted_equilibrated_flux_dirichlet( ...
                mesh,p,p,u_full,zero_full,I_eye(2,2), ...
                lambda_h,I_intval(0),flux_options);
        assert(isfinite(I_sup(rho_eig)) && I_inf(rho_eig) >= 0);
        assert(isfinite(I_sup( ...
            eig_flux_info.equilibrium_defect_norm_upper)));

        % A conservative mock lower bound for lambda_2 is sufficient for
        % this assembly-level certificate smoke test.
        L2_mock = I_intval(0.9*lambda_h_all(2));
        mu_upper = I_intval(I_sup(lambda_h));
        assert(I_inf(L2_mock-mu_upper) > 0);
        eta_eig = L2_mock/(L2_mock-mu_upper) ...
            *I_intval(I_sup(rho_eig));
        q2_eig = eta_eig^2/L2_mock;
        assert(I_sup(q2_eig) < 1);
    end

    data = struct();
    data.i = i;
    data.lambda_bounds = I_intval(lambda_h_all(:));
    data.lambda_h = lambda_h;
    data.lambda1_lower = I_intval(0.9*lambda_h_all(1));
    data.u_h = u;
    data.v_h = v;
    data.K = K;
    data.M = M;
    data.A_e = Ae;
    data.A_ee = Aee;
    data.P_e = Pe;
    data.P_ee = Pee;
    % Positive mock discretization errors: double mode is a sign/assembly
    % test, not a proof certificate.
    data.eps_a = 2e-2;
    data.eps_0 = 1e-3;
    data.flux_residual = eta_flux;

    [lower,upper,info] = residual_hessian_enclosure(data);
    assert(isfinite(I_inf(lower)) && isfinite(I_sup(upper)));
    assert(I_inf(lower) <= I_sup(upper));
    assert(I_inf(info.C_plus) > 0);
    if i == 1
        assert(I_sup(info.C_minus) == 0);
        assert(I_sup(info.upper_residual_correction) == 0);

        % The cell driver evaluates these quantities through transported
        % 2-by-2 gradient covariances.  Check that the scalar-form API is
        % algebraically identical to the matrix API at a point triangle.
        forms = struct();
        forms.mass_u = u'*M*u;
        forms.alpha_h = alpha_h;
        forms.c_h = (u'*Aee*u)/forms.mass_u;
        forms.ell_h = -u'*Ae*v+alpha_h*(u'*M*v);
        forms.b_h_vv = v'*K*v-lambda_h*(v'*M*v);
        forms.grad_u_sq = u'*K*u;
        forms.grad_v_sq = v'*K*v;
        forms.l2_v_sq = v'*M*v;
        scalar_data = data;
        scalar_data.scalar_forms = forms;
        [scalar_lower,scalar_upper,scalar_info] = ...
            residual_hessian_enclosure(scalar_data);
        assert(abs(I_inf(scalar_lower)-I_inf(lower)) < 1e-8);
        assert(abs(I_sup(scalar_upper)-I_sup(upper)) < 1e-8);
        assert(scalar_info.used_scalar_forms);

        % For i=1 the upper resolvent constant is C_minus=0, so the upper
        % Hessian enclosure is valid without constructing a shifted flux.
        upper_data = rmfield(scalar_data,'flux_residual');
        upper_data.upper_only = true;
        [upper_lower,upper_upper,upper_info] = ...
            residual_hessian_enclosure(upper_data);
        assert(isinf(I_inf(upper_lower)) && I_inf(upper_lower) < 0);
        assert(abs(I_sup(upper_upper)-I_sup(upper)) < 1e-8);
        assert(I_sup(upper_info.R) == 0);
    else
        assert(I_inf(info.C_minus) > 0);
        assert(I_inf(info.upper_residual_correction) > 0);
    end

    fprintf(['i=%d lambda_h=%.10g eta_flux=%.3e ', ...
        'ddlambda in [%.10g, %.10g], C-=%g C+=%g\n'], ...
        i,I_mid(lambda_h),I_sup(eta_flux),I_inf(lower),I_sup(upper), ...
        I_mid(info.C_minus),I_mid(info.C_plus));
end

fprintf('test_residual_hessian_estimator_double: PASS (INTERVAL_MODE=%d)\n', ...
    INTERVAL_MODE);
end


function A = shape_matrix(P,Axx,Axy,Ayy)
A = P(1,1)*Axx+(P(1,2)+P(2,1))*Axy+P(2,2)*Ayy;
end


function mesh = make_uniform_triangle_mesh(apex,n)
index = zeros(n+1,n+1);
nodes = zeros((n+1)*(n+2)/2,2);
cursor = 1;
for i = 0:n
    for j = 0:(n-i)
        index(i+1,j+1) = cursor;
        nodes(cursor,:) = (i/n)*[1,0]+(j/n)*apex;
        cursor = cursor+1;
    end
end

elements = zeros(n^2,3);
cursor = 1;
for i = 0:(n-1)
    for j = 0:(n-1-i)
        elements(cursor,:) = [index(i+1,j+1), ...
            index(i+2,j+1),index(i+1,j+2)];
        cursor = cursor+1;
        if i+j <= n-2
            elements(cursor,:) = [index(i+2,j+1), ...
                index(i+2,j+2),index(i+1,j+2)];
            cursor = cursor+1;
        end
    end
end
elements = elements(1:cursor-1,:);

all_edges = [elements(:,[1,2]);elements(:,[2,3]);elements(:,[3,1])];
all_edges = sort(all_edges,2);
[edges,~,edge_ids] = unique(all_edges,'rows');
counts = accumarray(edge_ids,1);
boundary_edges = edges(counts==1,:);

mesh = struct();
mesh.nodes = nodes;
mesh.elements = elements;
mesh.edges = edges;
mesh.boundary_edges = boundary_edges;
end
