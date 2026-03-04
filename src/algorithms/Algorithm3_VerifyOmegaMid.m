function [is_verified, results, diagnostics] = Algorithm3_VerifyOmegaMid( ...
    conjecture_type, cell_def_file, verbose, cell_range)
% ALGORITHM3_VERIFYOMEGAMID: Verify conjecture in Omega_mid region
%
% Inputs (all required):
%   conjecture_type: 'J1' or 'J2'
%   cell_def_file: path to cell definition CSV
%   verbose: logical
%   cell_range: [] or [start,end]   (indices in the cell_def_file order)
%
% Outputs:
%   is_verified, results, diagnostics
%
% IMPORTANT:
%   - No checkpoint files.
%   - No CSV log, no resume, no skipping.
%   - This function ONLY computes the requested range and returns results.

global INTERVAL_MODE;

% Initialize INTLAB mode if needed
if isempty(INTERVAL_MODE)
    INTERVAL_MODE = 0;
end

%% ========================================================================
%  STEP 2-1: Read cell definitions
%  ========================================================================

persistent CELLDEF_CACHE_KEY CELLDEF_CACHE_TABLE

% --- normalize path (string -> char) ---
fpath = char(cell_def_file);

% --- if relative path / on MATLAB path, resolve it ---
if ~isfile(fpath)
    alt = which(fpath);
    if ~isempty(alt)
        fpath = alt;
    end
end

% --- file stat for cache key ---
d = dir(fpath);
if isempty(d)
    error('Algorithm3_VerifyOmegaMid:CellDefNotFound', 'cell_def_file not found: %s', fpath);
end
d = d(1);

cache_key = sprintf('%s|%.15g|%d', fpath, d.datenum, d.bytes);

% --- cache load ---
if isempty(CELLDEF_CACHE_KEY) || ~strcmp(CELLDEF_CACHE_KEY, cache_key) || isempty(CELLDEF_CACHE_TABLE)
    CELLDEF_CACHE_TABLE = readtable(fpath);
    CELLDEF_CACHE_KEY   = cache_key;
end

cell_data = CELLDEF_CACHE_TABLE;
n_cells_total = height(cell_data);

% Determine cell range to process (indices in the file order)
if isempty(cell_range)
    cell_start = 1;
    cell_end   = n_cells_total;
else
    cell_start = max(1, floor(cell_range(1)));
    cell_end   = min(floor(cell_range(2)), n_cells_total);
end
if cell_end < cell_start
    tmp = cell_start; cell_start = cell_end; cell_end = tmp;
end
n_cells = cell_end - cell_start + 1;

%% ========================================================================
%  STEP 2-3: Process each cell using verify_J_positive
%  ========================================================================
n_verified = 0;
n_failed   = 0;

tic_total = tic;

% Allocate outputs for the requested range (in order)
all_results = struct();
all_results.cell_index   = (cell_start:cell_end).';
all_results.cell_id      = NaN(n_cells, 1);
all_results.verified     = false(n_cells, 1);
all_results.J_lower      = -Inf(n_cells, 1);
all_results.compute_time = NaN(n_cells, 1);
all_results.status       = strings(n_cells, 1);
all_results.note         = strings(n_cells, 1);


for i_cell = cell_start:cell_end
    local_k = i_cell - cell_start + 1;

    tic_cell = tic;
    current_cell = cell_data(i_cell, :);

    % Build the input struct expected by verify_J_positive
    cell_struct = struct();
    cell_struct.i = current_cell.i;
    cell_struct.x_inf = num2str(current_cell.x_inf, '%.17g');
    cell_struct.x_sup = num2str(current_cell.x_sup, '%.17g');
    cell_struct.theta_inf = num2str(current_cell.theta_inf, '%.17g');
    cell_struct.theta_sup = num2str(current_cell.theta_sup, '%.17g');
    cell_struct.mesh_size_upper = current_cell.mesh_size_upper;
    cell_struct.fem_order_upper = current_cell.fem_order_upper;
    cell_struct.mesh_size_lower_cr = current_cell.mesh_size_lower_cr;
    cell_struct.isLG = current_cell.isLG;

    if current_cell.isLG == 1
        cell_struct.mesh_size_lower_LG  = current_cell.mesh_size_lower_LG;
        cell_struct.fem_order_lower_LG  = current_cell.fem_order_lower_LG;
    else
        cell_struct.mesh_size_lower_LG = NaN;
        cell_struct.fem_order_lower_LG = NaN;
    end

    % Compute
    status = "ok";
    note   = "";
    try
        [is_cell_verified, J_lower, ~] = verify_J_positive(conjecture_type, cell_struct);
    catch ME
        if verbose
            fprintf('[Cell index %d] ERROR: %s\n', i_cell, ME.message);
            fprintf('%s\n', getReport(ME, 'extended', 'hyperlinks', 'on'));
        end
        J_lower = -Inf;
        is_cell_verified = false;
        status = "error";
        note   = string(ME.message);
    end

    cell_time = toc(tic_cell);

    % Store
    all_results.cell_id(local_k)      = double(current_cell.i);
    all_results.verified(local_k)     = logical(is_cell_verified);
    all_results.J_lower(local_k)      = J_lower;
    all_results.compute_time(local_k) = double(cell_time);
    all_results.status(local_k)       = status;
    all_results.note(local_k)         = note;

    if is_cell_verified
        n_verified = n_verified + 1;
    else
        n_failed = n_failed + 1;
        if verbose
            fprintf('[Cell index %d] FAILED: J_lower = %.6e <= 0\n', i_cell, J_lower);
        end
    end

    % Progress
    if verbose && mod(local_k, 100) == 0
        elapsed = toc(tic_total);
        rate = local_k / max(elapsed, eps);
        eta_seconds = (n_cells - local_k) / max(rate, eps);

        fprintf('[Progress] %d/%d (%.1f%%), Verified: %d, Failed: %d, ETA: %.1f min\n', ...
            local_k, n_cells, 100*local_k/n_cells, n_verified, n_failed, eta_seconds/60);
    end
end

time_elapsed = toc(tic_total);
is_verified = (n_failed == 0);

results = struct();
results.region          = 'Omega_mid';
results.conjecture_type = conjecture_type;
results.n_cells         = n_cells;
results.n_verified      = n_verified;
results.n_failed        = n_failed;
results.is_verified     = is_verified;
results.cell_results    = all_results;
results.time_elapsed    = time_elapsed;

if n_verified > 0
    verified_J = all_results.J_lower(all_results.verified);
    try
        results.J_min  = I_inf(min(verified_J));
        results.J_max  = I_sup(max(verified_J));
    catch
        results.J_min  = min(verified_J);
        results.J_max  = max(verified_J);
    end
    results.J_mean = mean(verified_J);
else
    results.J_min  = NaN;
    results.J_max  = NaN;
    results.J_mean = NaN;
end

diagnostics = struct();
diagnostics.cell_def_file = cell_def_file;
diagnostics.verbose       = verbose;
diagnostics.cell_range    = [cell_start, cell_end];
diagnostics.INTERVAL_MODE = INTERVAL_MODE;
diagnostics.time_elapsed  = time_elapsed;

end
