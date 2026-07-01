function rows = fem_lv_jxx_sweep(varargin)
%FEM_LV_JXX_SWEEP  Sweep FEM/Liu--Vejchodsky Jxx lower bounds at one triangle.
%
% Paper quantities evaluated here:
%   1. calc_ddlami_lower_bound evaluates the certified FEM
%      Dirichlet--Neumann estimate for ddot(lambda_1), i.e. the theorem
%      lower bound
%
%          D_{1,N}^{h,low}
%          - 2 * lambda_{N+1}/(lambda_{N+1}-lambda_1)
%              * (U_{1,M}^{h,up} - S_{1,N}^{h,low})
%          <= ddot(lambda_1).
%
%      The Dirichlet and Neumann eigenspace errors in this estimate are the
%      Liu--Vejchodsky cluster/subspace errors computed by
%      calc_grad_error_bounds.
%
%   2. This script then evaluates the paper's Omega_up x-convexity formula
%
%          Jxx_k = 0.5*y*ddot(lambda_1) + Rxx_k(x,y),   k = 1,2,
%
%      using the explicit Rxx_k terms from Algorithm2_VerifyOmegaUp.
%
% The default is an interval run.  Use 'interval_mode', false only for local
% smoke tests; non-interval mode cannot certify positivity.

p = inputParser();
addParameter(p, 'N_list', [1 3 6 10 17 30]);
addParameter(p, 'M_list', [2 5 7 10 20 50]);
addParameter(p, 'N_LG', 20);
addParameter(p, 'N_rho', 80);
addParameter(p, 'fem_ord_LG', 3);
addParameter(p, 'x', 0.5);
addParameter(p, 'y', NaN);
addParameter(p, 'interval_mode', true);
addParameter(p, 'output_dir', '');
addParameter(p, 'gmsh_command', '');
addParameter(p, 'intlab_path', '');
parse(p, varargin{:});
opt = p.Results;

project_root = fileparts(fileparts(mfilename('fullpath')));
if isempty(opt.output_dir)
    opt.output_dir = fullfile(project_root, 'results', 'fem_lv_jxx_sweep');
end
if ~exist(opt.output_dir, 'dir')
    mkdir(opt.output_dir);
end

configure_project(project_root, opt);

global INTERVAL_MODE
if INTERVAL_MODE
    rigorous_flag = 'interval_certified';
else
    rigorous_flag = 'double_smoke_test';
end

xI = scalar_to_interval(opt.x);
if isnan(opt.y)
    yI = sqrt(I_intval('3'))/2;
else
    yI = scalar_to_interval(opt.y);
end

base_triangle = I_intval([0, 0, 1, 0, xI, yI]);
target_triangle = base_triangle;
e_direction = [I_intval('1'), I_intval('0')];

run_stamp = datestr(now, 'yyyymmdd_HHMMSS');
csv_file = fullfile(opt.output_dir, ['fem_lv_jxx_sweep_' run_stamp '.csv']);
summary_file = fullfile(opt.output_dir, ['summary_' run_stamp '.md']);

header = ['timestamp,rigorous_flag,interval_mode,x_mid,y_mid,N,M,N_LG,N_rho,fem_ord_LG,ok,elapsed_sec,' ...
          'lambda1_inf,lambda1_sup,lambda1_width,dlambda_inf,dlambda_sup,' ...
          'ddlambda_lower,ddlambda_upper,ddlambda_width,J1_lower,J2_lower,Rxx_J1_inf,Rxx_J2_inf,' ...
          'D_lower,D_upper,S_lower,U_upper,tail_gap_upper,Q_tail_upper,G_upper,' ...
          'eta_phi_1,eta_phi_used_max,eta_psi_used_max,' ...
          'dir_clusters,neumann_clusters,min_dir_proj_margin,min_neu_proj_margin,error_message'];
write_text(csv_file, [header newline]);

rows = {};
row_count = 0;

fprintf('FEM/Liu--Vejchodsky Jxx sweep\n');
fprintf('  rigorous_flag=%s\n', rigorous_flag);
fprintf('  x=%.17g, y=%.17g\n', I_mid(xI), I_mid(yI));
fprintf('  N_LG=%d, N_rho=%d, fem_ord_LG=%d\n', opt.N_LG, opt.N_rho, opt.fem_ord_LG);
fprintf('  output=%s\n\n', csv_file);

for a = 1:numel(opt.N_list)
    N = opt.N_list(a);
    for b = 1:numel(opt.M_list)
        M = opt.M_list(b);
        fprintf('[N=%d, M=%d] ', N, M);
        t0 = tic();
        err_msg = '';
        ok = false;
        data = empty_row_data();

        try
            [lami, dlami, ddlambda_lower, diagnostics] = calc_ddlami_lower_bound( ...
                1, base_triangle, target_triangle, e_direction, ...
                N, opt.N_LG, opt.N_rho, opt.fem_ord_LG, M);

            Rxx1 = omega_up_Rxx('J1', xI, yI);
            Rxx2 = omega_up_Rxx('J2', xI, yI);
            Jxx1 = I_intval('0.5')*yI*ddlambda_lower + Rxx1;
            Jxx2 = I_intval('0.5')*yI*ddlambda_lower + Rxx2;

            dir_ids = diagnostics.dirichlet_cluster_ids;
            neu_ids = diagnostics.neumann_cluster_ids;

            data.lambda1_inf = lower(lami);
            data.lambda1_sup = upper(lami);
            data.lambda1_width = data.lambda1_sup - data.lambda1_inf;
            data.dlambda_inf = lower(dlami);
            data.dlambda_sup = upper(dlami);
            data.ddlambda_lower = lower(ddlambda_lower);
            data.ddlambda_upper = upper(diagnostics.ddlam_upper);
            data.ddlambda_width = data.ddlambda_upper - data.ddlambda_lower;
            data.J1_lower = lower(Jxx1);
            data.J2_lower = lower(Jxx2);
            data.Rxx_J1_inf = lower(Rxx1);
            data.Rxx_J2_inf = lower(Rxx2);
            data.D_lower = lower(diagnostics.D_lower);
            data.D_upper = upper(diagnostics.D_upper);
            data.S_lower = lower(diagnostics.S_lower);
            data.U_upper = upper(diagnostics.U_upper);
            data.tail_gap_upper = upper(diagnostics.tail_gap_upper);
            data.Q_tail_upper = upper(diagnostics.Q_tail_upper);
            data.G_upper = upper(diagnostics.G_upper);
            data.eta_phi_1 = upper(diagnostics.eta_phi_1);
            data.eta_phi_used_max = max_upper(diagnostics.eta_phi_cluster(dir_ids));
            data.eta_psi_used_max = max_upper(diagnostics.eta_psi_cluster(neu_ids));
            data.dir_clusters = encode_clusters(diagnostics.dirichlet_clusters, dir_ids);
            data.neumann_clusters = encode_clusters(diagnostics.mu_clusters, neu_ids);
            data.min_dir_proj_margin = min_margin(diagnostics.dir_proj_lower, diagnostics.dir_proj_radius);
            data.min_neu_proj_margin = min_margin(diagnostics.neu_proj_lower, diagnostics.neu_proj_radius);
            ok = true;

            fprintf('OK dd=%.6g J1=%.6g J2=%.6g S=%.6g U=%.6g etaD=%.3g etaN=%.3g\n', ...
                data.ddlambda_lower, data.J1_lower, data.J2_lower, ...
                data.S_lower, data.U_upper, data.eta_phi_used_max, data.eta_psi_used_max);
        catch ME
            err_msg = ME.message;
            fprintf('ERR %s\n', err_msg);
        end

        elapsed_sec = toc(t0);
        row_count = row_count + 1;
        rows{row_count, 1} = data; %#ok<AGROW>
        append_csv_row(csv_file, opt, rigorous_flag, xI, yI, N, M, ok, elapsed_sec, data, err_msg);
    end
end

write_summary(summary_file, csv_file, opt, rigorous_flag);
fprintf('\nWrote:\n  %s\n  %s\n', csv_file, summary_file);

end


function configure_project(project_root, opt)
global INTERVAL_MODE gmsh_command mesh_path

if isempty(opt.gmsh_command)
    env_gmsh = getenv('GMSH_COMMAND');
    if ~isempty(env_gmsh)
        gmsh_command = env_gmsh;
    elseif exist('/usr/bin/gmsh', 'file')
        gmsh_command = '/usr/bin/gmsh';
    elseif exist('/opt/homebrew/Cellar/gmsh/4.15.0/bin/gmsh', 'file')
        gmsh_command = '/opt/homebrew/Cellar/gmsh/4.15.0/bin/gmsh';
    else
        gmsh_command = 'gmsh';
    end
else
    gmsh_command = opt.gmsh_command;
end

mesh_path = fullfile(project_root, 'src', 'mesh');

if opt.interval_mode
    INTERVAL_MODE = 1;
    intlab_path = opt.intlab_path;
    if isempty(intlab_path)
        env_intlab = getenv('INTLAB_PATH');
        if ~isempty(env_intlab)
            intlab_path = env_intlab;
        else
            intlab_path = fullfile(project_root, 'Intlab_V12');
        end
    end
    addpath(intlab_path);
    evalc('startintlab');
else
    INTERVAL_MODE = 0;
end

addpath(fullfile(project_root, 'src', 'algorithms'));
addpath(fullfile(project_root, 'src', 'fem'));
addpath(fullfile(project_root, 'src', 'mesh'));
addpath(fullfile(project_root, 'src', 'interval'));
addpath(fullfile(project_root, 'src'));
addpath(fullfile(project_root, 'src', 'lib', 'VFEM2D', 'lib_eigenvalue_bound'));
addpath(fullfile(project_root, 'src', 'lib', 'VFEM2D_revised'));
addpath(fullfile(project_root, 'src', 'lib', 'veigs'));
addpath(fullfile(project_root, 'inputs'));
addpath(fullfile(project_root, 'results'));
addpath(fullfile(project_root, 'tests'));
addpath(fullfile(project_root, 'scripts_run'));
addpath(project_root);
end


function xI = scalar_to_interval(x)
if isa(x, 'intval')
    xI = x;
else
    xI = I_intval(sprintf('%.17g', x));
end
end


function RxxI = omega_up_Rxx(conjecture_type, xI, yI)
% Explicit Rxx_k in Lemma Ji-xx-yy-simpler-nopA, as used by Algorithm 2.
r1 = sqrt(xI.^2 + yI.^2);
r2 = sqrt((xI-1).^2 + yI.^2);
P  = 1 + r1 + r2;

switch upper(conjecture_type)
    case 'J1'
        RxxI = -(I_pi^2./(4*yI)) .* ( ...
            (xI./r1 + (xI-1)./r2).^2 + ...
            P.*(yI.^2).*(1./(r1.^3) + 1./(r2.^3)) );
    case 'J2'
        Cstar = 4*I_pi^2 / (3 + sqrt(I_pi*sqrt(I_intval('3'))))^2;
        Q = P + sqrt(2*I_pi*yI);
        RxxI = -(Cstar./yI) .* ( ...
            (xI./r1 + (xI-1)./r2).^2 + ...
            Q.*(yI.^2).*(1./(r1.^3) + 1./(r2.^3)) );
    otherwise
        error('Unknown conjecture_type: %s', conjecture_type);
end
end


function data = empty_row_data()
fields = {'lambda1_inf','lambda1_sup','lambda1_width','dlambda_inf','dlambda_sup', ...
          'ddlambda_lower','ddlambda_upper','ddlambda_width','J1_lower','J2_lower', ...
          'Rxx_J1_inf','Rxx_J2_inf','D_lower','D_upper','S_lower','U_upper', ...
          'tail_gap_upper','Q_tail_upper','G_upper','eta_phi_1','eta_phi_used_max', ...
          'eta_psi_used_max','min_dir_proj_margin','min_neu_proj_margin'};
for k = 1:numel(fields)
    data.(fields{k}) = NaN;
end
data.dir_clusters = '';
data.neumann_clusters = '';
end


function append_csv_row(csv_file, opt, rigorous_flag, xI, yI, N, M, ok, elapsed_sec, data, err_msg)
fid = fopen(csv_file, 'a');
if fid < 0
    error('Cannot open CSV file for append: %s', csv_file);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, ['%s,%s,%d,%.17g,%.17g,%d,%d,%d,%d,%d,%d,%.17g,' ...
              '%.17g,%.17g,%.17g,%.17g,%.17g,' ...
              '%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,' ...
              '%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,' ...
              '%.17g,%.17g,%.17g,%s,%s,%.17g,%.17g,%s\n'], ...
    csv_escape(datestr(now, 'yyyy-mm-dd HH:MM:SS')), csv_escape(rigorous_flag), logical(opt.interval_mode), ...
    I_mid(xI), I_mid(yI), N, M, opt.N_LG, opt.N_rho, opt.fem_ord_LG, ok, elapsed_sec, ...
    data.lambda1_inf, data.lambda1_sup, data.lambda1_width, data.dlambda_inf, data.dlambda_sup, ...
    data.ddlambda_lower, data.ddlambda_upper, data.ddlambda_width, data.J1_lower, data.J2_lower, ...
    data.Rxx_J1_inf, data.Rxx_J2_inf, data.D_lower, data.D_upper, data.S_lower, data.U_upper, ...
    data.tail_gap_upper, data.Q_tail_upper, data.G_upper, data.eta_phi_1, data.eta_phi_used_max, ...
    data.eta_psi_used_max, csv_escape(data.dir_clusters), csv_escape(data.neumann_clusters), ...
    data.min_dir_proj_margin, data.min_neu_proj_margin, csv_escape(err_msg));
end


function write_summary(summary_file, csv_file, opt, rigorous_flag)
fid = fopen(summary_file, 'w');
if fid < 0
    error('Cannot write summary file: %s', summary_file);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '# FEM/Liu--Vejchodsky Jxx Sweep\n\n');
fprintf(fid, 'This run evaluates the FEM Dirichlet--Neumann lower bound for `ddot(lambda_1)` and then `Jxx_k = 0.5*y*ddot(lambda_1) + Rxx_k` at one triangle.\n\n');
fprintf(fid, '- rigorous flag: `%s`\n', rigorous_flag);
fprintf(fid, '- N list: `%s`\n', mat2str(opt.N_list));
fprintf(fid, '- M list: `%s`\n', mat2str(opt.M_list));
fprintf(fid, '- mesh/order: `N_LG=%d`, `N_rho=%d`, `fem_ord_LG=%d`\n', opt.N_LG, opt.N_rho, opt.fem_ord_LG);
fprintf(fid, '- full CSV: `%s`\n\n', csv_file);
fprintf(fid, 'Use the CSV columns `J1_lower` and `J2_lower` to identify the first pair `(N,M)` certifying positivity.  The diagnostic columns `S_lower`, `U_upper`, `eta_phi_used_max`, `eta_psi_used_max`, and `min_dir_proj_margin` show whether the obstruction is the Dirichlet projection lower bound or the Neumann complement.\n');
end


function s = encode_clusters(clusters, ids)
parts = cell(1, numel(ids));
for k = 1:numel(ids)
    idx = clusters{ids(k)};
    parts{k} = sprintf('%d-%d', idx(1), idx(end));
end
s = strjoin(parts, ';');
end


function v = min_margin(norm_lower, radius)
if isempty(norm_lower)
    v = NaN;
    return;
end
vals = zeros(numel(norm_lower), 1);
for k = 1:numel(norm_lower)
    vals(k) = lower(norm_lower(k)) - upper(radius(k));
end
v = min(vals);
end


function v = max_upper(x)
if isempty(x)
    v = NaN;
else
    vals = zeros(numel(x), 1);
    for k = 1:numel(x)
        vals(k) = upper(x(k));
    end
    v = max(vals);
end
end


function y = lower(x)
y = I_inf(x);
if ~isreal(y), y = real(y); end
end


function y = upper(x)
y = I_sup(x);
if ~isreal(y), y = real(y); end
end


function write_text(path, text)
fid = fopen(path, 'w');
if fid < 0
    error('Cannot write file: %s', path);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', text);
end


function s = csv_escape(x)
if isstring(x), x = char(x); end
if islogical(x), x = char(string(x)); end
if isnumeric(x), x = num2str(x); end
if isempty(x), x = ''; end
x = strrep(char(x), '"', '""');
s = ['"' x '"'];
end
