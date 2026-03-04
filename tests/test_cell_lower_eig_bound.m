function test_cell_lower_eig_bound()
% Minimal smoke test for cell_lower_eig_bound (CR + LG)

% ---- add paths (adjust if you place this file elsewhere) ----
this = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(this, 'src')));   % repo has src/...
addpath(genpath(this));                   % in case cell_lower_eig_bound.m is elsewhere

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
cell_base.mesh_size_lower_cr = 0.08;

% only lambda_1 lower bound
cell_base.neig = 1;

fprintf('INTERVAL_MODE=%d\n', INTERVAL_MODE);

% ============================================================
% CR branch
% ============================================================
cell_cr = cell_base;
cell_cr.isLG = 0;

tic;
lb_cr = cell_lower_eig_bound(cell_cr);
t_cr = toc;

fprintf('CR : lb = %s   time = %.3f s\n', mat2str(lb_cr, 16), t_cr);
assert(all(isfinite(lb_cr)) && all(lb_cr > 0), 'CR lower bound invalid');

% ============================================================
% LG branch
% ============================================================
cell_lg = cell_base;
cell_lg.isLG = 1;
cell_lg.mesh_size_lower_LG  = 0.5;  % very coarse mesh (as you observed: tris=4)
cell_lg.fem_order_lower_LG  = 2;    % P2

tic;
lb_lg = cell_lower_eig_bound(cell_lg);
t_lg = toc;

fprintf('LG : lb = %s   time = %.3f s\n', mat2str(lb_lg, 16), t_lg);
assert(all(isfinite(lb_lg)) && all(lb_lg > 0), 'LG lower bound invalid');

fprintf('LG/CR ratio = %.6f\n', lb_lg(1) / lb_cr(1));
fprintf('OK\n');
end