function [is_verified, results, diagnostics] = Algorithm2_VerifyOmegaUp( ...
    conjecture_type, eps_up, Nx, Ny, Ny_axis, N_spec, mesh_params)
% ALGORITHM2_VERIFYOMEGAUP
% Verification in Ω_up (near the equilateral triangle) — Step 1 of the paper:
%
% Inputs (ALL REQUIRED):
%   conjecture_type : 'J1' or 'J2'
%   eps_up          : epsilon_up in Ω_up definition
%   Nx              : number of x-subintervals for R_ij covering
%   Ny              : number of y-subintervals for R_ij covering
%   Ny_axis         : number of y-subintervals for I_m^y on x=1/2
%   N_spec          : spectral truncation index N_spec
%   mesh_params     : struct for FEM/LG settings
%
% Outputs:
%   is_verified, results, diagnostics
%
% Date: 2026-01-15

global INTERVAL_MODE;
if isempty(INTERVAL_MODE)
    INTERVAL_MODE = 1;  % recommend interval mode for certification
end

% ----------------------------
% Ensure consistent types (as in the original parsing)
% ----------------------------
conjecture_type = char(conjecture_type);

% ----------------------------
% Basic constants / equilateral apex
% ----------------------------
p0 = [I_intval('0.5'), sqrt(I_intval('3'))/2];

eps_rect = eps_up;   % used for the R_ij cover around p0 

% Store diagnostics
diagnostics = struct();
diagnostics.INTERVAL_MODE = INTERVAL_MODE;
diagnostics.conjecture_type = conjecture_type;
diagnostics.eps_up = eps_up;
diagnostics.eps_rect = eps_rect;
diagnostics.Nx = Nx;
diagnostics.Ny = Ny;
diagnostics.Ny_axis = Ny_axis;
diagnostics.N_spec = N_spec;
diagnostics.mesh_params = mesh_params;

fprintf('================================================================\n');
fprintf('ALGORITHM 2 (Step 1): Ω_up verification (paper-matched algorithms)\n');
fprintf('================================================================\n');
fprintf('Conjecture: %s\n', conjecture_type);
fprintf('eps_up = %.17g\n', I_mid(eps_up));
fprintf('Nx=%d, Ny=%d (rect grid), Ny_axis=%d (axis intervals), N_spec=%d\n', Nx, Ny, Ny_axis, N_spec);
fprintf('Interval mode: %d\n\n', INTERVAL_MODE);

% =========================================================================
% STEP (1-2) Algorithm: certify sign of ∂²J_k/∂x² on Ω_up via rectangles R_ij
% =========================================================================
fprintf('--- STEP (1-2): Certify sign of ∂²J/∂x² over Ω_up (rectangular cover)\n');
dx_report = certify_sign_d2Jdx2_over_OmegaUp(conjecture_type, eps_rect, Nx, Ny, N_spec, mesh_params);
fprintf('    min lower bound over all R_ij: %.17g\n', I_inf(dx_report.min_lower_bound));
fprintf('    all cells certified positive: %s\n\n', mat2str(dx_report.all_positive));

% =========================================================================
% STEP (1-3) Algorithm: certify sign of ∂²J_k/∂y² on symmetry axis x=1/2
% =========================================================================
fprintf('--- STEP (1-3): Certify sign of ∂²J/∂y² on x=1/2 (axis subdivision)\n');
dy_report = certify_sign_d2Jdy2_on_axis(conjecture_type, eps_up, Ny_axis, N_spec, mesh_params);
fprintf('    min lower bound over all I_m^y: %.17g\n', I_inf(dy_report.min_lower_bound));
fprintf('    all intervals certified positive: %s\n\n', mat2str(dy_report.all_positive));

% =========================================================================
% Final decision
% =========================================================================
is_verified = dx_report.all_positive && dy_report.all_positive;

fprintf('================================================================\n');
fprintf('Ω_up verified: %s\n', mat2str(is_verified));
fprintf('================================================================\n\n');

% ----------------------------
% Pack results
% ----------------------------
results = struct();
results.region = 'Omega_up';
results.p0 = p0;

results.step_1_2 = dx_report;
results.step_1_3 = dy_report;

results.is_verified = is_verified;

end


% ========================================================================
% ALGORITHM (paper): Certification of the sign of ∂²J_k/∂x² (k=1,2)
% ========================================================================
function report = certify_sign_d2Jdx2_over_OmegaUp(conjecture_type, eps_rect, Nx, Ny, N_spec, mesh_params)
% Cover (a neighborhood of) Ω_up by rectangles R_ij using the paper grid:
%   p_ij = (1/2 + 2*eps_rect*i/Nx,  sqrt(3)/2 - eps_rect*j/Ny)
% Partition of [1/2, 1/2 + 2*eps_rect] x [sqrt(3)/2 - eps_rect, sqrt(3)/2],
% which covers Ω_up. Certify a rigorous *lower bound* on ∂²J/∂x² over each R_ij.

x0    = I_intval('0.5');
y_top = sqrt(I_intval('3'))/2;

x_nodes = x0    + 2*eps_rect*(0:Nx)/Nx;   % length Nx+1, spans [1/2, 1/2+2*eps]
y_nodes = y_top - eps_rect*(0:Ny)/Ny;     % length Ny+1, descending in [sqrt(3)/2-eps, sqrt(3)/2]

L  = I_zeros(Ny, Nx);                   % L(j,i) = lower bound on cell (i,j)
ok = true(Ny, Nx);

failures = {};                          % kept for compatibility

% -------------------------------------------------------------
% Per-cell CSV append (crash-safe in MATLAB):
%   open->write->close on every iteration (acts as flush)
% -------------------------------------------------------------
do_save    = true;
results_dir = 'results';

if isstruct(mesh_params)
    if isfield(mesh_params,'save_intermediate')
        do_save = logical(mesh_params.save_intermediate);
    end
    if isfield(mesh_params,'results_dir')
        results_dir = char(string(mesh_params.results_dir));
    end
end

csv_file = '';
if do_save
    if ~exist(results_dir,'dir'), mkdir(results_dir); end
    csv_file = fullfile(results_dir, sprintf('%s_OmegaUp_step1_2_cells.csv', conjecture_type));
    ensureCsvHeader(csv_file, ...
        'conjecture,step,i,j,x_lo_inf,x_hi_sup,y_lo_inf,y_hi_sup,L_lower,ok,min_so_far,run_timestamp');
end

min_so_far = +Inf;

for i = 1:Nx
    for j = 1:Ny
        % Cell R_ij: x in [x_i, x_{i+1}], y in [y_{j+1}, y_j]
        fprintf('cell verification:[i,j]=[%d,%d]\n', i,j);
        x_lo = x_nodes(i);
        x_hi = x_nodes(i+1);

        y_hi = y_nodes(j);
        y_lo = y_nodes(j+1);

        % Anchor point p_ij = (x_i, y_j) (top-left corner)
        p_anchor = [x_lo, y_hi];

        % Interval box for this rectangle
        box = struct('x', [x_lo, x_hi], 'y', [y_lo, y_hi]);

        % Main certification for this cell
        Lij = I_intval(NaN); % default: means "not computed / failed"

        if box_disjoint_from_Omega(box)
            fprintf('  -> skipped (disjoint from Omega)\n');
            Lij = I_intval(Inf);    % safe sentinel for skipped cells (won't reduce min)
        else
            Lij = lower_bound_d2Jdx2_on_rectangle(conjecture_type, p_anchor, box, N_spec, mesh_params);
        end

        % Store
        L(j,i)  = Lij;
        ok(j,i) = (Lij > 0);

        % Rigorous numeric lower bound to log
        Llower = scalarInf(Lij);
        if ~isnan(Llower)
            min_so_far = min(min_so_far, Llower);
        end

        % Append immediately (open->write->close)
        if do_save
            ts = datestr(now,'yyyy-mm-dd HH:MM:SS');
            appendCsvRow(csv_file, sprintf( ...
                '%s,%s,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17e,%d,%.17e,%s', ...
                conjecture_type, '1-2', i, j, ...
                scalarInf(x_lo), scalarSup(x_hi), ...
                scalarInf(y_lo), scalarSup(y_hi), ...
                Llower, double(Llower > 0), min_so_far, ts));
        end
    end
end

report = struct();
report.algorithm      = 'certify-d2Jdx2';
report.Nx             = Nx;
report.Ny             = Ny;
report.eps_rect       = eps_rect;
report.L              = L;
report.ok             = ok;
report.all_positive   = all(ok(:));
report.min_lower_bound = min(L(:));
report.failures       = failures;

% =============================================================
% local helpers
% =============================================================
function v = scalarInf(x)
    if isa(x,'intval'), v = I_inf(x); else, v = double(x); end
end

function v = scalarSup(x)
    if isa(x,'intval'), v = I_sup(x); else, v = double(x); end
end

function ensureCsvHeader(f, headerLine)
    fid = fopen(f,'a');
    if fid < 0
        warning('Cannot open CSV: %s', f);
        return;
    end

    fseek(fid, 0, 'eof');
    if ftell(fid) == 0
        fprintf(fid, '%s\n', headerLine);
    end
    fclose(fid);
end

function appendCsvRow(f, line)
    fid = fopen(f,'a');
    if fid < 0
        warning('Cannot append CSV: %s', f);
        return;
    end
    fprintf(fid, '%s\n', line);
    fclose(fid); % acts as "flush"
end

end

function Lij = lower_bound_d2Jdx2_on_rectangle(conjecture_type, p_anchor, box, N_spec, mesh_params)
% Implements the *single-cell* certification described in the paper algorithm.
%
% Output:
%   Lij: rigorous lower bound on ∂²J_k/∂x² over all p in this rectangle.

% Direction for ∂²/∂x² in the (x,y)-parameter space
e_direction = [I_intval('1'), I_intval('0')];

% Interval variables for the cell (x,y) ∈ box
xI = I_hull(box.x(1), box.x(2));
yI = I_hull(box.y(1), box.y(2));

% ------------------------------------------------------------------------
% [Lower bound of ∂²λ1/∂x²] using calc_ddlami_lower_bound at the anchor triangle
% ------------------------------------------------------------------------
base_triangle = I_intval([0, 0, 1, 0, box.x(1), box.y(1)]);
triangle      = I_intval([0, 0, 1, 0, xI, yI]);

[N_LG, N_rho, fem_ord_LG] = get_mesh_params_for_calc_ddlami(mesh_params);

[~, ~, ddlambda_lower] = calc_ddlami_lower_bound( ...
    1, base_triangle, triangle, e_direction, ...
    N_spec, N_LG, N_rho, fem_ord_LG);

% ------------------------------------------------------------------------
% [Lower bound of ∂²J_k/∂x²] via Lemma (Jkxx-simple):
%   ∂²J_k/∂x² = 0.5*y*∂²λ1/∂x² + R_xx^(k)(x,y)
% Take an interval enclosure and then I_inf() to get a rigorous lower bound.
% ------------------------------------------------------------------------
RxxI = interval_Rxx(conjecture_type, xI, yI);

JxxI = I_intval('0.5')*yI*ddlambda_lower + RxxI;
Lij = I_inf(JxxI);

end


% ========================================================================
% ALGORITHM (paper): Certification of the sign of ∂²J_k/∂y² on x=1/2
% ========================================================================
function report = certify_sign_d2Jdy2_on_axis(conjecture_type, eps_up, Ny_axis, N_spec, mesh_params)
% Axis interval:
%   I^y = { x=1/2, y ∈ [sqrt(3)/2 - eps_up, sqrt(3)/2] }
% Partition into Ny_axis subintervals and certify a lower bound on ∂²J/∂y²
% over each.

x_axis = I_intval('0.5');
y_top  = sqrt(I_intval('3'))/2;
y_min  = y_top - eps_up;

y_nodes = linspace(y_min, y_top, Ny_axis+1); % Ny_axis intervals

L  = I_zeros(Ny_axis, 1);
ok = true(Ny_axis, 1);

failures = {};

% -------------------------------------------------------------
% Per-interval CSV append (crash-safe in MATLAB):
%   open->write->close on every iteration
% -------------------------------------------------------------
do_save     = true;
results_dir = 'results';

if isstruct(mesh_params)
    if isfield(mesh_params,'save_intermediate')
        do_save = logical(mesh_params.save_intermediate);
    end
    if isfield(mesh_params,'results_dir')
        results_dir = char(string(mesh_params.results_dir));
    end
end

csv_file = '';
if do_save
    if ~exist(results_dir,'dir'), mkdir(results_dir); end
    csv_file = fullfile(results_dir, sprintf('%s_OmegaUp_step1_3_axis.csv', conjecture_type));
    ensureCsvHeader(csv_file, ...
        'conjecture,step,m,y_lo_inf,y_hi_sup,L_lower,ok,min_so_far,run_timestamp');
end

min_so_far = +Inf;

for m = Ny_axis:-1:1
    y_lo = y_nodes(m);
    y_hi = y_nodes(m+1);

    % Anchor point p_m = (1/2, y_m) (lower endpoint)
    p_anchor = [x_axis, y_lo];

    % Interval "cell" I_m^y: x fixed, y interval
    box = struct('x', x_axis, 'y', [y_lo, y_hi]);

    Lm = lower_bound_d2Jdy2_on_axis_interval(conjecture_type, p_anchor, box, N_spec, mesh_params);
    L(m)  = Lm;
    ok(m) = (Lm > 0);

    Llower = scalarInf(Lm);
    if ~isnan(Llower)
        min_so_far = min(min_so_far, Llower);
    end

    if do_save
        ts = datestr(now,'yyyy-mm-dd HH:MM:SS');
        appendCsvRow(csv_file, sprintf( ...
            '%s,%s,%d,%.17g,%.17g,%.17e,%d,%.17e,%s', ...
            conjecture_type, '1-3', m, ...
            scalarInf(y_lo), scalarSup(y_hi), ...
            Llower, double(Llower > 0), min_so_far, ts));
    end
end

report = struct();
report.algorithm       = 'certify-d2Jdy2';
report.Ny_axis         = Ny_axis;
report.eps_up          = eps_up;
report.L               = L;
report.ok              = ok;
report.all_positive    = all(ok);
report.min_lower_bound = min(L);
report.failures        = failures;

% =============================================================
% local helpers
% =============================================================
function v = scalarInf(x)
    if isa(x,'intval'), v = I_inf(x); else, v = double(x); end
end

function v = scalarSup(x)
    if isa(x,'intval'), v = I_sup(x); else, v = double(x); end
end

function ensureCsvHeader(f, headerLine)
    if exist(f,'file'), return; end
    fid = fopen(f,'w');
    if fid < 0
        warning('Cannot create CSV: %s', f);
        return;
    end
    fprintf(fid, '%s\n', headerLine);
    fclose(fid);
end

function appendCsvRow(f, line)
    fid = fopen(f,'a');
    if fid < 0
        warning('Cannot append CSV: %s', f);
        return;
    end
    fprintf(fid, '%s\n', line);
    fclose(fid); % acts as "flush"
end

end



function Lm = lower_bound_d2Jdy2_on_axis_interval(conjecture_type, p_anchor, box, N_spec, mesh_params)
% Implements the *single-interval* certification described in the paper algorithm.
%
% IMPORTANT CHANGE (requested):
%   Use calc_ddlami_lower_bound (instead of theorem3p1_lower_bound_ddlambda / lower_bound_dotlambda).
%
% Uses Lemma (Jkyy-simple):
%   ∂²J_k/∂y² = 0.5*y*∂²λ1/∂y² + ∂λ1/∂y + R_yy^(k).

% Direction for ∂/∂y and ∂²/∂y² in the (x,y)-parameter space
e_direction = [I_intval('0'), I_intval('1')];

xI = box.x; % constant 1/2, but keep in interface
yI = I_hull(box.y(1), box.y(2));

% ------------------------------------------------------------------------
% [Lower bound of ∂²λ1/∂y²] and [lower bound of ∂λ1/∂y]
% using calc_ddlami_lower_bound at the anchor triangle
% ------------------------------------------------------------------------
base_triangle = I_intval([0, 0, 1, 0, p_anchor(1), p_anchor(2)]);
triangle      = [0, 0, I_intval('1'), 0, xI, yI];

[N_LG, N_rho, fem_ord_LG] = get_mesh_params_for_calc_ddlami(mesh_params);

[~, dlambdaI, ddlambda_lower] = calc_ddlami_lower_bound( ...
    1, base_triangle, triangle, e_direction, ...
    N_spec, N_LG, N_rho, fem_ord_LG);

% Certified lower bound for ∂λ1/∂y over the anchor computation
dotlambda_lower = I_intval((I_inf(dlambdaI)));

% ------------------------------------------------------------------------
% [Lower bound of ∂²J_k/∂y²] via Lemma (Jkyy-simple)
% ------------------------------------------------------------------------
RyyI = interval_Ryy(conjecture_type, xI, I_intval(I_inf(yI)));

JyyI = I_intval('0.5')*yI*ddlambda_lower + dotlambda_lower + RyyI;
Lm = I_inf(JyyI);

end


% ========================================================================
% Explicit remainder terms R_xx^(k), R_yy^(k) from Lemma (Ji-xx-yy-simpler-nopA)
% ========================================================================
function RxxI = interval_Rxx(conjecture_type, xI, yI)
r1 = sqrt(xI.^2 + yI.^2);
r2 = sqrt((xI-1).^2 + yI.^2);
P  = 1 + r1 + r2;

switch upper(conjecture_type)
    case 'J1'
        RxxI = -(I_pi^2./(4*yI)) .* ( ...
            (xI./r1 + (xI-1)./r2).^2 + ...
            P.*(yI.^2).*(1./(r1.^3) + 1./(r2.^3)) );
    case 'J2'
        Cstar = 4*I_pi^2 / (3 + sqrt(I_pi*sqrt(3)))^2;
        Q = P + sqrt(2*I_pi*yI); % = 1+r1+r2 + sqrt(2πy)
        RxxI = -(Cstar./yI) .* ( ...
            (xI./r1 + (xI-1)./r2).^2 + ...
            Q.*(yI.^2).*(1./(r1.^3) + 1./(r2.^3)) );
    otherwise
        error('Unknown conjecture_type: %s (expected J1 or J2).', conjecture_type);
end
end


function RyyI = interval_Ryy(conjecture_type, xI, yI)
r1 = sqrt(xI.^2 + yI.^2);
r2 = sqrt((xI-1).^2 + yI.^2);
P  = 1 + r1 + r2;

switch upper(conjecture_type)
    case 'J1'
        RyyI = -(I_pi^2./(4*yI.^3)) .* ( ...
            (yI.^2).*P.*(xI.^2./(r1.^3) + (xI-1).^2./(r2.^3)) + ...
            (1 + xI.^2./r1 + (xI-1).^2./r2).^2 );
    case 'J2'
        Cstar = 4*I_pi^2 / (3 + sqrt(I_pi*sqrt(3)))^2;

        Q    = P + sqrt(2*I_pi*yI);
        Qy   = yI./r1 + yI./r2 + sqrt(I_pi./(2*yI));
        Qyy  = xI.^2./(r1.^3) + (xI-1).^2./(r2.^3) - 0.5*sqrt(I_pi/2).*yI.^(-3/2);

        RyyI = -(Cstar/2) .* ( ...
            (2./yI).*(Qy.^2 + Q.*Qyy) ...
            - (4./(yI.^2)).*(Q.*Qy) ...
            + (2./(yI.^3)).*(Q.^2) );
    otherwise
        error('Unknown conjecture_type: %s (expected J1 or J2).', conjecture_type);
end
end


% ========================================================================
% dotP / ddotP interval matrices (depend only on y and direction e=(a,b))
% ========================================================================
function dotP = interval_dotP(yI, e)
a = e(1); b = e(2);
dotP = [ I(0,0),      -a./yI ; ...
        -a./yI,  (-2*b)./yI ];
end

function ddotP = interval_ddotP(yI, e)
a = e(1); b = e(2);
ddotP = [ (2*a^2)./(yI.^2), (4*a*b)./(yI.^2) ; ...
          (4*a*b)./(yI.^2), (6*b^2)./(yI.^2) ];
end


% ========================================================================
% Mesh parameter extraction for calc_ddlami_lower_bound
% ========================================================================
function [N_LG, N_rho, fem_ord_LG] = get_mesh_params_for_calc_ddlami(mesh_params)
% Minimal helper to keep the Algorithm2 interface unchanged.
%
% Required fields for calc_ddlami_lower_bound:
%   N_LG, N_rho, fem_ord_LG
%
% If your project uses different field names, adapt them here (only here).

if isfield(mesh_params, 'N_LG')
    N_LG = mesh_params.N_LG;
elseif isfield(mesh_params, 'N_LG_points')
    N_LG = mesh_params.N_LG_points;
else
    error('mesh_params must contain N_LG (or N_LG_points).');
end

if isfield(mesh_params, 'N_rho')
    N_rho = mesh_params.N_rho;
else
    error('mesh_params must contain N_rho.');
end

if isfield(mesh_params, 'fem_ord_LG')
    fem_ord_LG = mesh_params.fem_ord_LG;
elseif isfield(mesh_params, 'fem_ord_lg')
    fem_ord_LG = mesh_params.fem_ord_lg;
else
    error('mesh_params must contain fem_ord_LG (or fem_ord_lg).');
end
end

function tf = box_disjoint_from_Omega(box)
% Return true only if the rectangle box is PROVABLY disjoint from
% Omega = {x^2+y^2 <= 1, x >= 1/2, y > 0}.

xI = I_hull(box.x(1), box.x(2));
yI = I_hull(box.y(1), box.y(2));


if I_sup(xI) < 0.5
    tf = true; return;
end

if I_sup(yI) <= 0
    tf = true; return;
end

sI = xI.^2 + yI.^2;
if I_inf(sI) > 1
    tf = true; return;
end

tf = false;
end