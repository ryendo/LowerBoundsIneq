function test_cell_lower_eig_bound()
% Minimal smoke test for cell_lower_eig_bound (CR + LG)

project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(project_root,'scripts_run'));
omega_up_all_prepare_worker(project_root,'interval');

% ---- choose interval / non-interval mode ----
global INTERVAL_MODE
INTERVAL_MODE = 1;   % true: veigs (interval), false: eigs (faster)

% ---- one representative cell (your example) ----
cell_base = struct();

cell_base.x_inf  = '0.5';
cell_base.x_sup  = '0.5';
cell_base.theta_inf = I_pi/3;
cell_base.theta_sup = I_pi/3;

% For speed in a test you can try 0.06~0.10
cell_base.mesh_size_lower_cr = 0.02;

% only lambda_1 lower bound
cell_base.neig = 1;

fprintf('INTERVAL_MODE=%d\n', INTERVAL_MODE);

% ============================================================
% CR branch
% ============================================================
cell_cr = cell_base;
cell_cr.isLG = 0;

tic;
[lb_cr,diag_cr] = cell_lower_eig_bound(cell_cr);
t_cr = toc;

fprintf('CR : lb = %s   time = %.3f s\n', mat2str(lb_cr, 16), t_cr);
assert(all(isfinite(lb_cr)) && all(lb_cr > 0), 'CR lower bound invalid');
assert(diag_cr.conditions_verified && strcmp(diag_cr.method,'CR'));
lb_cr_legacy = local_legacy_cr_lambda1(cell_cr);
assert(abs(lb_cr-lb_cr_legacy) ...
    <= 1e-10*max(1,abs(lb_cr_legacy)), ...
    'Two-value CR solve changed the legacy four-value lambda_1 bound.');

% ============================================================
% LG branch
% ============================================================
cell_lg = cell_base;
cell_lg.isLG = 1;
cell_lg.mesh_size_lower_LG  = 0.1;  % very coarse mesh (as you observed: tris=4)
cell_lg.fem_order_lower_LG  = 2;

tic;
[lb_lg,diag_lg] = cell_lower_eig_bound(cell_lg);
t_lg = toc;

fprintf('LG : lb = %s   time = %.3f s\n', mat2str(lb_lg, 16), t_lg);
assert(all(isfinite(lb_lg)) && all(lb_lg > 0), 'LG lower bound invalid');
assert(diag_lg.conditions_verified && strcmp(diag_lg.method,'LG'));
assert(I_inf(diag_lg.separation_margin) > 0);
assert(I_sup(diag_lg.mu_upper) < 0);
assert(strcmp(diag_lg.shift_method,'CR-Liu'));
assert(~diag_lg.fallback_used);
assert(diag_lg.timing.cr_shift_seconds > 0);
assert(diag_lg.timing.lg_setup_seconds >= 0);

% One of the two legacy rows that had been serialized as +Inf must now
% pass the complete strict-LG path with finite bounds for both functionals.
T = readtable(fullfile(project_root,'inputs','cell_def.csv'));
row = T(find(T.i==138748,1),:);
cell_mid = struct( ...
    'i',row.i, ...
    'x_inf',num2str(row.x_inf,'%.17g'), ...
    'x_sup',num2str(row.x_sup,'%.17g'), ...
    'theta_inf',num2str(row.theta_inf,'%.17g'), ...
    'theta_sup',num2str(row.theta_sup,'%.17g'), ...
    'mesh_size_lower_cr',row.mesh_size_lower_cr, ...
    'isLG',row.isLG, ...
    'mesh_size_lower_LG',row.mesh_size_lower_LG, ...
    'fem_order_lower_LG',row.fem_order_lower_LG);
[both_ok,J_both,both_diag] = verify_J_both_positive(cell_mid);
assert(both_ok);
assert(isfinite(J_both.J1) && J_both.J1 > 0);
assert(isfinite(J_both.J2) && J_both.J2 > 0);
assert(I_inf(both_diag.eigenvalue.separation_margin) > 0);
assert(I_sup(both_diag.eigenvalue.mu_upper) < 0);
assert(~both_diag.eigenvalue.fallback_used);
assert(strcmp(both_diag.eigenvalue.shift_method,'CR-Liu'));
assert(isempty(both_diag.eigenvalue.fallback_reason));
assert(both_diag.eigenvalue.timing.cr_shift_seconds > 0);

% A second legacy-Inf row must use the same CR--Liu LG shift policy.
row_fast = T(find(T.i==152387,1),:);
cell_fast = struct( ...
    'i',row_fast.i, ...
    'x_inf',num2str(row_fast.x_inf,'%.17g'), ...
    'x_sup',num2str(row_fast.x_sup,'%.17g'), ...
    'theta_inf',num2str(row_fast.theta_inf,'%.17g'), ...
    'theta_sup',num2str(row_fast.theta_sup,'%.17g'), ...
    'mesh_size_lower_cr',row_fast.mesh_size_lower_cr, ...
    'isLG',row_fast.isLG, ...
    'mesh_size_lower_LG',row_fast.mesh_size_lower_LG, ...
    'fem_order_lower_LG',row_fast.fem_order_lower_LG);
[fast_ok,J_fast,fast_diag] = verify_J_both_positive(cell_fast);
assert(fast_ok);
assert(isfinite(J_fast.J1) && J_fast.J1 > 0);
assert(isfinite(J_fast.J2) && J_fast.J2 > 0);
assert(~fast_diag.eigenvalue.fallback_used);
assert(strcmp(fast_diag.eigenvalue.shift_method,'CR-Liu'));
assert(fast_diag.eigenvalue.timing.cr_shift_seconds > 0);

fprintf('LG/CR ratio = %.6f\n', lb_lg(1) / lb_cr(1));
fprintf('OK\n');
end


function lower = local_legacy_cr_lambda1(cell_data)
x = I_intval(cell_data.x_sup);
theta = I_intval(cell_data.theta_sup);
a = x;
b = x*tan(theta);
mesh = make_mesh_by_gmsh(a,b,cell_data.mesh_size_lower_cr);
is_boundary = ismember(mesh.edges,mesh.boundary_edges,'rows');
tri_by_edge = find_tri2edge(mesh.elements,mesh.edges);
[Mcr,Kcr] = create_matrix_crouzeix_raviart( ...
    mesh.elements,mesh.edges,mesh.nodes,tri_by_edge);
dofs = 1:size(mesh.edges,1);
dofs(find(is_boundary > 0)) = []; %#ok<FNDSB>
Mcr = Mcr(dofs,dofs);
Kcr = Kcr(dofs,dofs);
[values,indices] = veigs(Kcr,Mcr,4,'sm');
location = find(indices(:) == 1);
assert(numel(location) == 1);
lambda1_discrete = values(location);
Ch = I_intval('0.1893') ...
    *find_mesh_hmax(mesh.nodes,mesh.edges);
lower = I_inf(lambda1_discrete/(1+lambda1_discrete*Ch^2));
end
