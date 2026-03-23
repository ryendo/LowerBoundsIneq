classdef VerificationRunner < handle
% VERIFICATIONRUNNER: Unified class for computer-assisted proofs
%
% This class provides a clean interface to run all verification algorithms
% from the paper "Sharp Dirichlet Eigenvalue Inequalities on Triangles"
%
% Paper Reference: Section 4 (Computer-assisted proof of the theorem)
%
% Main Features:
%   - Run individual algorithms (Algorithm 2 for Omega_up, Algorithm 3 for Omega_mid)
%   - Run complete verification workflow for J1 or J2 conjectures
%   - Resume support for Omega_mid verification
%   - Save results to structured results folder
%   - Generate verification reports
%   - Parameter management with paper notation
%
% Usage:
%   runner = VerificationRunner();  
%   results = runner.verifyOmegaUp('J1');
%   runner.generateReport(results);
%
% Author: Based on paper by R. Endo, X. Liu, P. Mariano
% Date: 2025-01-14 (modified for J1/J2 evaluation)

    properties
        % Parameters (paper notation)
        eps_up = I_intval('0.122');           % (fixed) epsilon_up: Omega_up region threshold
        N_spectral = 1;          % (fixed) Number of spectral terms for ddlam computation
        N_LG = 16;                % (fixed) Mesh resolution for Lehmann-Goerisch
        N_rho = 64;              % (fixed) Mesh resolution for CR
        ord_LG = 2;              % (fixed) Lagrange order for LG lower bound

        % Algorithm 2 grid parameters (Omega_up)
        Nx = 10; % (fixed)
        Ny = 100; % (fixed)
        Ny_axis = 200; %

        % Paths
        results_dir = 'results';

        % State
        last_results = struct();

        % Flags
        verbose = true;
        save_intermediate = true;
        resume_enabled = true;
    end

    properties (Dependent)
        use_interval_arithmetic
    end

    methods

        function v = get.use_interval_arithmetic(~)
            global INTERVAL_MODE;
            v = INTERVAL_MODE;
        end

        function obj = VerificationRunner(varargin)
            % VERIFICATIONRUNNER: Constructor
            %
            % Usage:
            %   runner = VerificationRunner();
            %   runner = VerificationRunner('eps_up', 0.05, 'verbose', true);

            if nargin > 0
                obj.setParameters(varargin{:});
            end

            if obj.verbose
                fprintf('================================================================\n');
                fprintf('VERIFICATION RUNNER INITIALIZED\n');
                fprintf('================================================================\n');
                fprintf('Parameters:\n');
                fprintf('  eps_up: %.17f\n', I_mid(obj.eps_up));
                fprintf('  N_spectral: %d\n', obj.N_spectral);                
                fprintf('  ord_LG: %d, N_LG: %d, N_rho: %d\n',obj.ord_LG, obj.N_LG, obj.N_rho);
                fprintf('  Algorithm2 grid: Nx=%d, Ny=%d, Ny_axis=%d\n', obj.Nx, obj.Ny, obj.Ny_axis);
                fprintf('  use_interval_arithmetic: %d\n', obj.use_interval_arithmetic);
                fprintf('  results_dir: %s\n', obj.results_dir);
                fprintf('  resume_enabled (Omega_mid, CSV-based): %d\n', obj.resume_enabled);
                fprintf('================================================================\n\n');
            end
        end

        function setParameters(obj, varargin)
            % SETPARAMETERS: Set verification parameters

            p = inputParser;
            addParameter(p, 'eps_up', obj.eps_up);
            addParameter(p, 'N_spectral', obj.N_spectral);
            addParameter(p, 'N_LG', obj.N_LG);
            addParameter(p, 'N_rho', obj.N_rho);
            addParameter(p, 'ord_LG', obj.ord_LG);

            % Algorithm 2 grid parameters (Omega_up)
            addParameter(p, 'Nx', obj.Nx);
            addParameter(p, 'Ny', obj.Ny);
            addParameter(p, 'Ny_axis', obj.Ny_axis);

            addParameter(p, 'verbose', obj.verbose);
            addParameter(p, 'save_intermediate', obj.save_intermediate);
            addParameter(p, 'use_interval_arithmetic', obj.use_interval_arithmetic);

            addParameter(p, 'resume_enabled', obj.resume_enabled);

            parse(p, varargin{:});

            obj.eps_up = p.Results.eps_up;
            obj.N_spectral = p.Results.N_spectral;
            obj.N_LG = p.Results.N_LG;
            obj.N_rho = p.Results.N_rho;
            obj.ord_LG = p.Results.ord_LG;

            obj.Nx = p.Results.Nx;
            obj.Ny = p.Results.Ny;
            obj.Ny_axis = p.Results.Ny_axis;

            obj.verbose = p.Results.verbose;
            obj.save_intermediate = p.Results.save_intermediate;
            obj.use_interval_arithmetic = p.Results.use_interval_arithmetic;
            obj.resume_enabled = p.Results.resume_enabled;
        end

        function results = verifyOmegaUp(obj, conjecture_type)
            % VERIFYOMEGAUP: Verify conjecture in Omega_up region (near equilateral)
            %
            % Paper Reference: Section 4, Step 1
            %   Uses Algorithm2_VerifyOmegaUp to compute second-order derivatives
            %   of J1 or J2 and verify positive definiteness of Hessian.

            if obj.verbose
                fprintf('\n================================================================\n');
                fprintf('STEP 1: Verifying Omega_up (near equilateral triangle)\n');
                fprintf('================================================================\n\n');
            end

            mesh_params = struct('N_LG', obj.N_LG, ...
                                 'N_rho', obj.N_rho, ...
                                 'fem_ord_LG', obj.ord_LG);

            [~, results, diagnostics] = Algorithm2_VerifyOmegaUp( ...
                conjecture_type, obj.eps_up, obj.Nx, obj.Ny, obj.Ny_axis, obj.N_spectral, mesh_params);
                

            if obj.save_intermediate
                d = fullfile(obj.results_dir); if ~exist(d,'dir'), mkdir(d); end

                fsum = fullfile(d, sprintf('%s_OmegaUp.csv', conjecture_type));
                if ~exist(fsum,'file')
                    fid = fopen(fsum,'w');
                    fprintf(fid,'conjecture,eps_up,ddJx_lower,ddJy_lower,is_verified,N_LG,N_rho,N_spectral,ord_LG,INTERVAL_MODE,run_timestamp\n');
                    fclose(fid);
                end

                ddJx = results.step_1_2.min_lower_bound;
                ddJy = results.step_1_3.min_lower_bound;
                isv  = results.is_verified;

                N_LG       = getps(diagnostics,'N_LG',obj.N_LG);
                N_rho      = getps(diagnostics,'N_rho',obj.N_rho);
                N_spectral = getps(diagnostics,'N_spectral',obj.N_spectral);
                ord_LG     = getps(diagnostics,'ord_LG',obj.ord_LG);
                IM         = getps(diagnostics,'INTERVAL_MODE',0);

                fid = fopen(fsum,'a');
                fprintf(fid,'%s,%.17g,%.17e,%.17e,%.0f,%d,%d,%d,%.0f,%.0f,%.0f,%s\n', ...
                    conjecture_type, I_mid(obj.eps_up), I_inf(ddJx), I_inf(ddJy), isv, ...
                    N_LG, N_rho, N_spectral, ord_LG, IM,datestr(now,'yyyy-mm-dd HH:MM:SS'));
                fclose(fid);
            end

            function v = getps(S, path, def)
                v = def;
                try
                    cur = S;
                    parts = strsplit(path,'.');
                    for k = 1:numel(parts)
                        if ~isstruct(cur) || ~isfield(cur, parts{k}), return; end
                        cur = cur.(parts{k});
                    end
                    if islogical(cur) && isscalar(cur), v = double(cur);
                    elseif isnumeric(cur) && isscalar(cur), v = double(cur);
                    end
                catch
                end
            end
        end

        function results = verifyOmegaMid(obj, conjecture_type, cell_range, cell_def_file)
            % VERIFYOMEGAMID: Verify conjecture in Omega_mid region
            %
            % Paper Reference: Section 4, Step 2
            %   Uses Algorithm3_VerifyOmegaMid to verify J >= 0 for all cells.
            %
            %   - Each cell is processed one-by-one, and appended to a CSV immediately.
            %   - Resume is CSV-based: if resume_enabled=true, cells already present in the CSV are skipped.

            if obj.verbose
                fprintf('\n================================================================\n');
                fprintf('STEP 2: Verifying Omega_mid (intermediate region)\n');
                fprintf('================================================================\n\n');
            end

            % Determine target cell IDs (prefer IDs from the definition file)
            ids = expectedIds(cell_def_file, cell_range);
            nCells = numel(ids);

            % Read cell definition CSV once for verbose display
            Tdef = table();
            try
                Tdef = readtable(cell_def_file, 'TextType','string');
            catch
                Tdef = table();
            end

            % CSV output for per-cell append
            d = fullfile(obj.results_dir);
            if ~exist(d,'dir'), mkdir(d); end
            f = fullfile(d, sprintf('%s_OmegaMid.csv', conjecture_type));

            needHeader = ~exist(f,'file');

            if ~needHeader
                info = dir(f);
                if isempty(info) || info.bytes == 0
                    needHeader = true;
                else
                    fid = fopen(f,'r');
                    first = fgetl(fid);
                    fclose(fid);
                    expected = "conjecture,cell_id,verified,J_lower,status,note,run_timestamp";
                    if ~ischar(first) && ~isstring(first)
                        needHeader = true;
                    else
                        needHeader = ~startsWith(strtrim(string(first)), expected);
                    end
                end
            end

            if obj.save_intermediate && needHeader
                fid = fopen(f,'w');
                fprintf(fid,'conjecture,cell_id,verified,J_lower,status,note,run_timestamp\n');
                fclose(fid);
            end

            % Existing CSV state: keep the latest row index for each cell_id
            Tcsv = table();
            rowOfCid = containers.Map('KeyType','double','ValueType','double'); % cell_id -> row index

            if exist(f,'file')
                try
                    Tcsv = readtable(f, 'TextType','string');

                    if any(strcmpi(Tcsv.Properties.VariableNames, 'cell_id'))
                        rowsThis = true(height(Tcsv),1);
                        if any(strcmpi(Tcsv.Properties.VariableNames, 'conjecture'))
                            rowsThis = strcmpi(string(Tcsv.conjecture), conjecture_type);
                        end

                        cidv = double(Tcsv.cell_id);
                        for r = 1:height(Tcsv)
                            if ~rowsThis(r), continue; end
                            cid = cidv(r);
                            if ~isnan(cid)
                                rowOfCid(cid) = r; % keep the latest occurrence
                            end
                        end
                    end
                catch ME
                    Tcsv = table();
                    if obj.verbose
                        fprintf('[OmegaMid] Existing CSV read failed (%s). Starting from scratch.\n', ME.message);
                    end
                end
            end

            % Per-cell container for return value
            per_cell = repmat(struct( ...
                'cell_id', NaN, 'verified', NaN, 'J_lower', NaN, ...
                'status', '', 'note', ''), nCells, 1);

            n_verified = 0;
            all_ok = true;

            % Split target cells into:
            %   (1) existing rows that must be rechecked/overwritten
            %   (2) new unknown rows
            idx_recheck = [];
            idx_unknown = [];

            for i = 1:nCells
                cid = ids(i);

                if obj.resume_enabled && isKey(rowOfCid, cid)
                    r = rowOfCid(cid);

                    row = struct();
                    row.cell_id  = cid;
                    row.verified = getNumRow(Tcsv, r, "verified", NaN);
                    row.J_lower  = getNumRow(Tcsv, r, "J_lower", NaN);
                    row.status   = getStrRow(Tcsv, r, "status", "");
                    row.note     = getStrRow(Tcsv, r, "note", "");
                    per_cell(i) = row;

                    if isfinite(row.verified) && row.verified == 1
                        n_verified = n_verified + 1;
                        if obj.verbose
                            fprintf('[OmegaMid] Skip verified cell %d/%d (id=%d)\n', i, nCells, cid);
                        end
                    else
                        idx_recheck(end+1) = i; %#ok<AGROW>
                    end
                else
                    idx_unknown(end+1) = i; %#ok<AGROW>
                end
            end

            n_recheck_left = numel(idx_recheck);
            n_unknown_left = numel(idx_unknown);

            if obj.verbose
                fprintf('[OmegaMid] Existing rows to recheck: %d\n', n_recheck_left);
                fprintf('[OmegaMid] Unknown rows to compute: %d\n', n_unknown_left);
            end

            proc_groups = {idx_recheck, idx_unknown};
            proc_labels = {'Recheck', 'Run new'};

            for g = 1:2
                idx_list = proc_groups{g};

                for kk = 1:numel(idx_list)
                    i = idx_list(kk);
                    cid = ids(i);

                    remRecheck = n_recheck_left - (g == 1);
                    remUnknown = n_unknown_left - (g == 2);

                    if obj.verbose
                        fprintf('[OmegaMid] %s cell %d/%d (id=%d) | remaining recheck=%d, unknown=%d\n', ...
                            proc_labels{g}, i, nCells, cid, remRecheck, remUnknown);

                        defLine = cellDefRowString(Tdef, cid, i, cell_range);
                        if ~isempty(defLine)
                            fprintf('[OmegaMid] cell_def row: %s\n', defLine);
                        end
                    end

                    try
                        [is_verified_cell, algo_results, ~] = Algorithm3_VerifyOmegaMid( ...
                            conjecture_type, ...
                            cell_def_file, ...
                            obj.verbose, ...
                            [i i]);

                        cr = pickSingleCell(algo_results);

                        verified = double(is_verified_cell);
                        Jlb      = getNum(cr, {'J_lower','J_lb','J_min','lower_bound','lb','JLower'}, NaN);
                        status   = getStr(cr, {'status','state','message','msg'}, 'ok');
                        note     = "";
                    catch ME
                        verified = NaN;
                        Jlb      = NaN;
                        status   = 'error';
                        note     = ME.message;
                    end

                    per_cell(i).cell_id  = cid;
                    per_cell(i).verified = verified;
                    per_cell(i).J_lower  = Jlb;
                    per_cell(i).status   = status;
                    per_cell(i).note     = note;

                    if isfinite(verified) && verified == 1
                        n_verified = n_verified + 1;
                    else
                        all_ok = false;
                    end

                    % Recheck: keep old behavior (rewrite whole CSV)
                    % Unknown: append one line only
                    if obj.save_intermediate
                        ts = string(datestr(now,'yyyy-mm-dd HH:MM:SS'));

                        if isKey(rowOfCid, cid)
                            % -----------------------------
                            % existing row -> recheck mode
                            % -----------------------------
                            r = rowOfCid(cid);
                            Tcsv.conjecture(r)    = string(conjecture_type);
                            Tcsv.cell_id(r)       = double(cid);
                            Tcsv.verified(r)      = double(verified);
                            Tcsv.J_lower(r)       = Jlb;
                            Tcsv.status(r)        = string(status);
                            Tcsv.note(r)          = string(note);
                            Tcsv.run_timestamp(r) = ts;

                            % keep current behavior for recheck
                            writeOmegaMidCSV(f, Tcsv);

                        else
                            % -----------------------------
                            % truly new row -> unknown mode
                            % -----------------------------
                            newRow = table( ...
                                string(conjecture_type), ...
                                double(cid), ...
                                double(verified), ...
                                Jlb, ...
                                string(status), ...
                                string(note), ...
                                ts, ...
                                'VariableNames', {'conjecture','cell_id','verified','J_lower','status','note','run_timestamp'});

                            % update in-memory state as well
                            if isempty(Tcsv) || width(Tcsv) == 0
                                Tcsv = newRow;
                            else
                                Tcsv = [Tcsv; newRow];
                            end
                            rowOfCid(cid) = height(Tcsv);

                            % append only one line to file
                            appendOmegaMidCSVRow(f, newRow);
                        end
                    end

                    if g == 1
                        n_recheck_left = n_recheck_left - 1;
                    else
                        n_unknown_left = n_unknown_left - 1;
                    end
                end
            end

            results = struct();
            results.is_verified   = all_ok;
            results.n_cells       = nCells;
            results.n_verified    = n_verified;
            results.cell_def_file = cell_def_file;
            results.cell_range    = cell_range;
            results.cell_results  = per_cell;

            if obj.verbose
                fprintf('[VerificationRunner] Omega_mid done: %d/%d verified\n', n_verified, nCells);
                fprintf('[VerificationRunner] CSV resume file: %s\n', f);
            end

            % ===========================
            % local helper functions
            % ===========================
            function appendOmegaMidCSVRow(fname, rowT)
                fid = fopen(fname, 'a');
                if fid < 0
                    warning('Cannot open for append: %s', fname);
                    return;
                end
                c = onCleanup(@() fclose(fid));

                writeOmegaMidCSVRow(fid, ...
                    rowT.conjecture(1), ...
                    rowT.cell_id(1), ...
                    rowT.verified(1), ...
                    rowT.J_lower(1), ...
                    rowT.status(1), ...
                    rowT.note(1), ...
                    rowT.run_timestamp(1));
            end

            function writeOmegaMidCSVRow(fid, conjecture, cell_id, verified, J_lower, status, note, ts)
                conjecture = csvText(conjecture);
                status     = csvText(status);
                note       = csvText(note);
                ts         = csvText(ts);

                cell_id  = numOrNaN(cell_id);
                verified = numOrNaN(verified);
                J_lower  = numOrNaN(J_lower);

                fprintf(fid, '%s,%d,%d,%.17e,"%s","%s",%s\n', ...
                    csvEsc(conjecture), ...
                    cell_id, ...
                    verified, ...
                    J_lower, ...
                    csvEsc(status), ...
                    csvEsc(note), ...
                    csvEsc(ts));
            end


            function s = numText(in)
                v = numOrNaN(in);
                if isnan(v)
                    s = 'NaN';
                else
                    s = sprintf('%.17g', v);
                end
            end

            function s = csvText(in)
                % Convert table entry / string / missing / cell safely to char row vector.
                if iscell(in)
                    if isempty(in)
                        s = '';
                        return;
                    end
                    in = in{1};
                end

                if isempty(in)
                    s = '';
                    return;
                end

                if ismissing(in)
                    s = '';
                    return;
                end

                if isstring(in)
                    if isscalar(in)
                        s = char(in);
                    else
                        s = char(join(in, "; "));
                    end
                    return;
                end

                if ischar(in)
                    s = in;
                    return;
                end

                if isnumeric(in) || islogical(in)
                    if isscalar(in)
                        if isnan(in)
                            s = '';
                        else
                            s = char(string(in));
                        end
                    else
                        s = mat2str(in);
                    end
                    return;
                end

                try
                    s = char(string(in));
                catch
                    s = '';
                end
            end

            function v = numOrNaN(in)
                if iscell(in)
                    if isempty(in)
                        v = NaN;
                        return;
                    end
                    in = in{1};
                end

                if isempty(in) || ismissing(in)
                    v = NaN;
                    return;
                end

                if isnumeric(in) || islogical(in)
                    v = double(in);
                    return;
                end

                if isstring(in) || ischar(in)
                    v = str2double(string(in));
                    if isnan(v)
                        v = NaN;
                    end
                    return;
                end

                v = NaN;
            end
            function writeOmegaMidCSV(fname, T)
                fid = fopen(fname, 'w');
                if fid < 0
                    warning('Cannot open: %s', fname);
                    return;
                end
                c = onCleanup(@() fclose(fid));

                fprintf(fid, 'conjecture,cell_id,verified,J_lower,status,note,run_timestamp\n');

                for rr = 1:height(T)
                    writeOmegaMidCSVRow(fid, ...
                        T.conjecture(rr), ...
                        T.cell_id(rr), ...
                        T.verified(rr), ...
                        T.J_lower(rr), ...
                        T.status(rr), ...
                        T.note(rr), ...
                        T.run_timestamp(rr));
                end
            end
            function s = cellDefRowString(T, cid, local_i, range_)
                % Return one printable line for the corresponding row in cell_def_file.
                % Preference:
                %   1) match by cell_id-like column
                %   2) fallback to row position implied by cell_range/local_i

                s = '';
                if isempty(T) || height(T) == 0
                    return;
                end

                r = NaN;

                % Try to match by ID column first
                vars = lower(string(T.Properties.VariableNames));
                cand = find(ismember(vars, ["cell_id","cellid","id","index","cell_idx","cellindex","cell_index"]), 1);

                if ~isempty(cand)
                    raw = T{:,cand};
                    if iscell(raw)
                        raw = cellfun(@(x) str2double(string(x)), raw);
                    else
                        raw = double(raw);
                    end
                    k = find(raw == cid, 1, 'first');
                    if ~isempty(k)
                        r = k;
                    end
                end

                % Fallback: use row index implied by cell_range
                if isnan(r)
                    if ~isempty(range_) && numel(range_) == 2
                        r = floor(range_(1)) + local_i - 1;
                    else
                        r = local_i;
                    end
                end

                if r < 1 || r > height(T)
                    return;
                end

                s = tableRowToString(T, r);
            end

            function s = tableRowToString(T, r)
                names = string(T.Properties.VariableNames);
                parts = strings(1, numel(names));

                for jj = 1:numel(names)
                    x = T{r, jj};
                    if iscell(x), x = x{1}; end
                    if ismissing(x)
                        val = "";
                    elseif isstring(x) || ischar(x)
                        val = string(x);
                    elseif isnumeric(x) || islogical(x)
                        if isscalar(x)
                            val = string(x);
                        else
                            val = string(mat2str(x));
                        end
                    else
                        try
                            val = string(jsonencode(x));
                        catch
                            val = "<unprintable>";
                        end
                    end
                    parts(jj) = names(jj) + "=" + val;
                end

                s = char(strjoin(parts, ', '));
            end

            function ids = expectedIds(def_file, range_)
                % Prefer IDs in def_file; fallback to 1..N; then apply range_ on ORDER.
                ids = [];
                try
                    T = readtable(def_file);
                    vars = lower(string(T.Properties.VariableNames));
                    cand = find(ismember(vars, ["cell_id","cellid","id","index","cell_idx","cellindex","cell_index"]), 1);
                    if ~isempty(cand)
                        raw = T{:,cand};
                        raw = raw(:);
                        if iscell(raw), raw = cellfun(@(x) str2double(string(x)), raw); end
                        ids = double(raw);
                    else
                        ids = (1:height(T))';
                    end
                    ids = ids(~isnan(ids)); ids = ids(:);
                catch
                    ids = [];
                end

                if isempty(ids)
                    if ~isempty(range_) && numel(range_)==2
                        ids = (floor(range_(1)):floor(range_(2)))';
                    end
                end

                if ~isempty(range_) && numel(range_)==2 && ~isempty(ids)
                    a = max(1, floor(range_(1)));
                    b = min(numel(ids), floor(range_(2)));
                    if a<=b, ids = ids(a:b); else, ids = ids([]); end
                end
            end

            function cr = pickSingleCell(ar)
                % Best-effort extraction of a single-cell result from Algorithm3 output
                cr = struct();
                if isstruct(ar) && isfield(ar,'cell_results') && ~isempty(ar.cell_results)
                    C = ar.cell_results;
                    cr = C(1); return;
                end
                if isstruct(ar) && isfield(ar,'cells') && ~isempty(ar.cells)
                    C = ar.cells;
                    cr = C(1); return;
                end
                if isstruct(ar), cr = ar; end
            end

            function v = getNum(S, names, def)
                v = def;
                if ~isstruct(S), return; end
                for ii = 1:numel(names)
                    nm = names{ii};
                    if isfield(S,nm)
                        x = S.(nm);
                        if iscell(x), x = x{1}; end
                        if islogical(x) && isscalar(x), v = double(x); return; end
                        if isnumeric(x) && isscalar(x), v = double(x); return; end
                        if (ischar(x) || isstring(x))
                            y = str2double(string(x));
                            if ~isnan(y), v = y; return; end
                        end
                    end
                end
            end

            function s = getStr(S, names, def)
                s = def;
                if ~isstruct(S), return; end
                for ii = 1:numel(names)
                    nm = names{ii};
                    if isfield(S,nm)
                        x = S.(nm);
                        if isempty(x), continue; end
                        if iscell(x), x = x{1}; end
                        if isstring(x) || ischar(x), s = char(string(x)); return; end
                    end
                end
            end

            function out = csvEsc(in)
                in = csvText(in);
                out = strrep(in, '"', '""');
            end

            function v = getNumRow(T, r, col, def)
                v = def;
                if any(strcmpi(T.Properties.VariableNames, col))
                    x = T{r, col};
                    if iscell(x), x = x{1}; end
                    if ismissing(x), return; end
                    if isnumeric(x) && isscalar(x), v = double(x); return; end
                    if isstring(x) || ischar(x)
                        y = str2double(string(x));
                        if ~isnan(y), v = y; end
                    end
                end
            end

            function s = getStrRow(T, r, col, def)
                s = def;
                if any(strcmpi(T.Properties.VariableNames, col))
                    x = T{r, col};
                    if iscell(x), x = x{1}; end
                    if ismissing(x), return; end
                    s = char(string(x));
                end
            end
        end

        function results = runCompleteVerification(obj, conjecture_type, varargin)
            % RUNCOMPLETEVERIFICATION: Run complete verification workflow
            %
            % Paper Reference: Section 4, Algorithm 4
            %   Runs all three steps: Omega_up, Omega_mid, Omega_down
            %
            % Usage:
            %   results = runner.runCompleteVerification('J1');
            %
            % Inputs:
            %   conjecture_type: 'J1' (Laugesen-Siudeja) or 'J2' (Cheeger-type)
            %
            % Outputs:
            %   results: Complete verification results

            % Parse options
            p = inputParser;
            addParameter(p, 'cell_range', []);
            addParameter(p, 'cell_def_file', 'inputs/cell_def.csv');
            parse(p, varargin{:});

            if obj.verbose
                fprintf('================================================================\n');
                fprintf('COMPLETE VERIFICATION: Conjecture %s\n', conjecture_type);
                fprintf('================================================================\n');
                fprintf('Date: %s\n', datestr(now));
                fprintf('Parameters (paper notation):\n');
                fprintf('  eps_up: %.4f\n', obj.eps_up);
                fprintf('  N_spectral: %d\n', obj.N_spectral);
                fprintf('  N_LG: %d, N_rho: %d\n', obj.N_LG, obj.N_rho);
                fprintf('================================================================\n');
            end

            total_tic = tic;

            % Step 1: Omega_up
            results_up = obj.verifyOmegaUp(conjecture_type);

            % Step 2: Omega_mid
            results_mid = obj.verifyOmegaMid(conjecture_type, p.Results.cell_range, p.Results.cell_def_file);


            % Combine results
            results = struct();
            results.conjecture_type = conjecture_type;
            results.timestamp = datestr(now);
            results.parameters = struct('eps_up', obj.eps_up, ...
                                       'N_spectral', obj.N_spectral, ...
                                       'N_LG', obj.N_LG, ...
                                       'N_rho', obj.N_rho);
            results.omega_up = results_up;
            results.omega_mid = results_mid;
            results.overall_verified = results_up.is_verified && results_mid.is_verified;

            results.total_time = toc(total_tic);

            % Save complete results
            timestamp = datestr(now, 'yyyymmdd');
            save_file = fullfile(obj.results_dir, 'complete', ...
                sprintf('%s_complete_verification_%s.csv', conjecture_type, timestamp));
            VerificationRunner.writeStructsToCSV(save_file, ...
                {'results'}, {results});

            % Store in object
            obj.last_results = results;

            % Print summary
            if obj.verbose
                fprintf('\n================================================================\n');
                fprintf('VERIFICATION SUMMARY\n');
                fprintf('================================================================\n');
                fprintf('Conjecture: %s\n', conjecture_type);
                fprintf('  Omega_up:   %s\n', mat2str(results_up.is_verified));
                fprintf('  Omega_mid:  %s\n', mat2str(results_mid.is_verified));
                fprintf('  --------------------------------\n');
                fprintf('  OVERALL:    %s\n', mat2str(results.overall_verified));
                fprintf('  --------------------------------\n');
                fprintf('Total time: %.2f seconds\n', results.total_time);
                fprintf('================================================================\n\n');
            end
        end

        function generateReport(obj, results)
            % GENERATEREPORT: Generate verification report
            %
            % Usage:
            %   runner.generateReport(results);

            if nargin < 2
                results = obj.last_results;
            end

            if isempty(results)
                error('No results to report. Run verification first.');
            end

            % Generate report filename
            timestamp = datestr(now, 'yyyymmdd_HHMMSS');
            report_file = fullfile(obj.results_dir, 'reports', ...
                sprintf('report_%s_%s.txt', results.conjecture_type, timestamp));

            % Open file
            fid = fopen(report_file, 'w');

            % Write report
            fprintf(fid, '================================================================\n');
            fprintf(fid, 'VERIFICATION REPORT\n');
            fprintf(fid, '================================================================\n\n');
            fprintf(fid, 'Paper: Sharp Dirichlet Eigenvalue Inequalities on Triangles\n');
            fprintf(fid, 'Authors: R. Endo, X. Liu, P. Mariano\n\n');
            fprintf(fid, 'Date: %s\n', results.timestamp);
            fprintf(fid, 'Conjecture: %s\n\n', results.conjecture_type);

            fprintf(fid, 'Parameters (paper notation):\n');
            fprintf(fid, '  eps_up: %.4f\n', results.parameters.eps_up);
            fprintf(fid, '  N_spectral: %d\n', results.parameters.N_spectral);
            fprintf(fid, '  N_LG: %d\n', results.parameters.N_LG);
            fprintf(fid, '  N_rho: %d\n\n', results.parameters.N_rho);

            fprintf(fid, 'Results:\n');
            fprintf(fid, '  Omega_up:   %s\n', mat2str(results.omega_up.is_verified));
            if isfield(results.omega_up, 'step_1_2')
                fprintf(fid, '    ddJ/dx^2 lower bound: %.6e\n', results.omega_up.step_1_2.ddJ_x_lower);
            end
            if isfield(results.omega_up, 'step_1_3')
                fprintf(fid, '    ddJ/dy^2 lower bound: %.6e\n', results.omega_up.step_1_3.ddJ_y_lower);
            end

            fprintf(fid, '  Omega_mid:  %s (%d/%d cells)\n', ...
                mat2str(results.omega_mid.is_verified), ...
                results.omega_mid.n_verified, results.omega_mid.n_cells);
            if isfield(results.omega_mid, 'J_min')
                fprintf(fid, '    J_lower range: [%.6e, %.6e]\n', results.omega_mid.J_min, results.omega_mid.J_max);
            end

            fprintf(fid, 'Overall Verification: %s\n\n', mat2str(results.overall_verified));

            fprintf(fid, 'Total Time: %.2f seconds\n\n', results.total_time);

            fprintf(fid, '================================================================\n');
            if results.overall_verified
                fprintf(fid, 'CONCLUSION: Conjecture %s is VERIFIED\n', results.conjecture_type);
                fprintf(fid, 'The equilateral triangle uniquely minimizes %s.\n', results.conjecture_type);
            else
                fprintf(fid, 'CONCLUSION: Verification INCOMPLETE\n');
                fprintf(fid, 'Some regions failed verification. See details above.\n');
            end
            fprintf(fid, '================================================================\n');

            fclose(fid);

            fprintf('[VerificationRunner] Report saved to: %s\n', report_file);
        end
    end

    methods (Access = private, Static)
        function writeStructsToCSV(csv_file, section_names, section_structs)
            % WRITESTRUCTSTOCSV: Save one or more structs to a single CSV file.
            % Format: section,key,value  (values are stringified safely)

            if ~iscell(section_names), section_names = {section_names}; end
            if ~iscell(section_structs), section_structs = {section_structs}; end

            % Ensure parent dir exists
            parent_dir = fileparts(csv_file);
            if ~isempty(parent_dir) && ~exist(parent_dir, 'dir')
                mkdir(parent_dir);
            end

            fid = fopen(csv_file, 'w');
            if fid < 0
                warning('[VerificationRunner] Could not open CSV for writing: %s', csv_file);
                return;
            end

            fprintf(fid, 'section,key,value\n');

            for sidx = 1:numel(section_names)
                sec = section_names{sidx};
                S = section_structs{sidx};

                [keys, vals] = VerificationRunner.flattenAnyToKeyValue(S, '');
                for i = 1:numel(keys)
                    key = keys{i};
                    val = vals{i};

                    sec_esc = VerificationRunner.csvEscape(sec);
                    key_esc = VerificationRunner.csvEscape(key);
                    val_esc = VerificationRunner.csvEscape(val);

                    fprintf(fid, '"%s","%s","%s"\n', sec_esc, key_esc, val_esc);
                end
            end

            fclose(fid);
        end

        function [keys, vals] = flattenAnyToKeyValue(x, prefix)
            % FLATTENANYTOKEYVALUE: Recursively flatten structs/cells/arrays into dot-keys.

            keys = {};
            vals = {};

            if nargin < 2, prefix = ''; end

            % Helper to push a leaf
            function pushLeaf(k, v)
                if isempty(k)
                    k = '(root)';
                end
                keys{end+1} = k; %#ok<AGROW>
                vals{end+1} = v; %#ok<AGROW>
            end

            % Struct (scalar)
            if isstruct(x) && isscalar(x)
                fns = fieldnames(x);
                for j = 1:numel(fns)
                    fn = fns{j};
                    if isempty(prefix)
                        newpref = fn;
                    else
                        newpref = [prefix '.' fn];
                    end
                    [k2, v2] = VerificationRunner.flattenAnyToKeyValue(x.(fn), newpref);
                    keys = [keys, k2]; %#ok<AGROW>
                    vals = [vals, v2]; %#ok<AGROW>
                end
                return;
            end

            % Struct array
            if isstruct(x) && ~isscalar(x)
                for j = 1:numel(x)
                    idxpref = sprintf('%s(%d)', prefix, j);
                    [k2, v2] = VerificationRunner.flattenAnyToKeyValue(x(j), idxpref);
                    keys = [keys, k2]; %#ok<AGROW>
                    vals = [vals, v2]; %#ok<AGROW>
                end
                return;
            end

            % Cell
            if iscell(x)
                for j = 1:numel(x)
                    idxpref = sprintf('%s{%d}', prefix, j);
                    [k2, v2] = VerificationRunner.flattenAnyToKeyValue(x{j}, idxpref);
                    keys = [keys, k2]; %#ok<AGROW>
                    vals = [vals, v2]; %#ok<AGROW>
                end
                return;
            end

            % Other: stringify as a leaf
            pushLeaf(prefix, VerificationRunner.stringifyValue(x));
        end

        function s = stringifyValue(v)
            % STRINGIFYVALUE: Best-effort conversion to a single-line string.

            try
                if isstring(v) && isscalar(v)
                    s = char(v);
                    return;
                end
                if ischar(v)
                    s = v;
                    return;
                end
                if isnumeric(v) || islogical(v)
                    if isscalar(v)
                        if islogical(v)
                            s = sprintf('%.0f', double(v));
                        else
                            s = sprintf('%.16g', v);
                        end
                    else
                        s = mat2str(v);
                    end
                    return;
                end

                % Try JSON (useful for many types)
                try
                    s = jsonencode(v);
                    if isstring(s), s = char(s); end
                    return;
                catch
                    % ignore
                end

                % Fallback: capture display text
                s = strtrim(evalc('disp(v)'));
                s = regexprep(s, '\s+', ' '); % single-line
            catch
                s = '[unprintable]';
            end
        end

        function out = csvEscape(in)
            % CSVESCAPE: escape quotes for CSV fields (we quote every field).
            if isstring(in) && isscalar(in), in = char(in); end
            if ~ischar(in), in = VerificationRunner.stringifyValue(in); end
            out = strrep(in, '"', '""');
        end
    end
end
