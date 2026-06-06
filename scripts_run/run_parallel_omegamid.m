function run_parallel_omegamid(conjecture_type, cell_range, nworkers, cell_def_file, out_dir)
%RUN_PARALLEL_OMEGAMID Parpool-driven verification of Omega_mid.
%   Splits cell_range into nworkers chunks; each worker writes to its own
%   subdirectory out_dir/worker_NN/. Resume is handled per-chunk by
%   VerificationRunner's CSV-based skip logic. After all chunks complete,
%   the merged result CSV is written to out_dir/<conj>_OmegaMid.csv.
%
% Inputs:
%   conjecture_type : 'J1' or 'J2'
%   cell_range      : [lo hi] (row indices in cell_def CSV). [] = all rows.
%   nworkers        : number of local parpool workers
%   cell_def_file   : path to inputs/cell_def_*.csv
%   out_dir         : parent results directory (workers live under it)

    if nargin < 5 || isempty(out_dir),       out_dir       = 'results_parallel'; end
    if nargin < 4 || isempty(cell_def_file), cell_def_file = 'inputs/cell_def_v12.csv'; end
    if nargin < 3 || isempty(nworkers),      nworkers      = 20; end

    % Client init: add paths, load INTLAB .mat (pre-generated). No delete.
    addpath(fileparts(mfilename("fullpath")));
    my_intlab_config_worker();

    % Determine total cell count from the input CSV.
    T = readtable(cell_def_file);
    n_rows = height(T);
    if nargin < 2 || isempty(cell_range)
        cell_range = [1, n_rows];
    end
    lo0 = cell_range(1); hi0 = min(cell_range(2), n_rows);
    total = hi0 - lo0 + 1;
    chunk = ceil(total / nworkers);

    if ~exist(out_dir,'dir'), mkdir(out_dir); end

    % Start / reuse parpool.
    p = gcp('nocreate');
    if isempty(p) || p.NumWorkers ~= nworkers
        if ~isempty(p), delete(p); end
        p = parpool('local', nworkers);
    end

    % Initialize INTLAB on every worker ONCE (before any parfor work).
    fprintf('[parallel] initializing INTLAB on %d workers...\n', nworkers);
    addpath(fullfile(fileparts(mfilename('fullpath'))));  % so workers see my_intlab_config_worker
    f_init = parfevalOnAll(gcp, @my_intlab_config_worker, 0);
    wait(f_init);  % CRITICAL: wait for per-worker INTLAB init before parfor, else some iterations fail to dispatch
    fprintf('[parallel] worker init done.\n');

    % Build per-worker chunk boundaries (precomputed to avoid closure pitfalls).
    los = zeros(1,nworkers);
    his = zeros(1,nworkers);
    for w = 1:nworkers
        los(w) = lo0 + (w-1)*chunk;
        his(w) = min(lo0 + w*chunk - 1, hi0);
    end

    fprintf('[parallel] %s on cells [%d, %d], %d workers, chunk=%d\n', ...
        conjecture_type, lo0, hi0, nworkers, chunk);
    t_start = tic;

    parfor w = 1:nworkers
        lo = los(w); hi = his(w);
        if lo > hi, continue; end
        wdir = fullfile(out_dir, sprintf('worker_%02d', w));
        if ~exist(wdir,'dir'), mkdir(wdir); end

        % Construct with no args to avoid setter on dependent use_interval_arithmetic.
        runner_w = VerificationRunner();
        runner_w.verbose           = false;
        runner_w.resume_enabled    = true;
        runner_w.save_intermediate = true;
        runner_w.results_dir       = wdir;

        fprintf('[worker %02d] %s cells [%d, %d] ->  %s\n', w, conjecture_type, lo, hi, wdir);
        try
            runner_w.verifyOmegaMid(conjecture_type, [lo hi], cell_def_file);
        catch ME
            fprintf(2, '[worker %02d] ERROR: %s\n', w, ME.message);
        end
    end

    fprintf('[parallel] parfor done in %.1f min\n', toc(t_start)/60);

    % Merge worker chunks into out_dir/<conj>_OmegaMid.csv
    merge_worker_chunks(out_dir, conjecture_type, nworkers);
end

function merge_worker_chunks(out_dir, conjecture_type, nworkers)
    fname = sprintf('%s_OmegaMid.csv', conjecture_type);
    merged_path = fullfile(out_dir, fname);
    fid_out = fopen(merged_path, 'w');
    fprintf(fid_out, 'conjecture,cell_id,verified,J_lower,status,note,run_timestamp\n');
    total = 0;
    for w = 1:nworkers
        chunk_path = fullfile(out_dir, sprintf('worker_%02d', w), fname);
        if ~exist(chunk_path, 'file'), continue; end
        fid_in = fopen(chunk_path, 'r');
        fgetl(fid_in);  % skip header
        while true
            line = fgetl(fid_in);
            if ~ischar(line), break; end
            if ~isempty(strtrim(line))
                fprintf(fid_out, '%s\n', line);
                total = total + 1;
            end
        end
        fclose(fid_in);
    end
    fclose(fid_out);
    fprintf('[merge] wrote %d rows -> %s\n', total, merged_path);
end
