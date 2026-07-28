classdef VerifyTriangleInequalities < handle
%VERIFYTRIANGLEINEQUALITIES  Main class for computer-assisted verification of
%   sharp Dirichlet Laplacian eigenvalue inequalities on planar triangles.
%
% Usage:
%   (1) Validate committed results without recomputing:
%           v = VerifyTriangleInequalities();
%           v.run();
%   (2) Reproduce the full computation from scratch:
%           run('scripts_run/run_omega_mid_unified_full.m')
%
% .compute(20) remains as a compatibility entry point, but delegates to the
% same strict unified wrapper.  That wrapper shares one CR/Liu or
% Liu--Lehmann--Goerisch spectral computation between J1 and J2, resumes
% validated chunk checkpoints, publishes the two CSVs transactionally, and
% finally calls .run().
%
% Paper: R. Endo, X. Liu, P. Mariano, "Sharp Dirichlet Eigenvalue
% Inequalities on Triangles".

    properties
        project_root   char
        inputs_dir     char
        results_dir    char
        scripts_dir    char
        raw_dir        char
        cell_def_file  char
        cell_def       table
        J1             table
        J2             table
    end

    methods
        function obj = VerifyTriangleInequalities(project_root)
            if nargin < 1 || isempty(project_root)
                project_root = fileparts(mfilename('fullpath'));
            end
            obj.project_root  = project_root;
            obj.inputs_dir    = fullfile(project_root, 'inputs');
            obj.results_dir   = fullfile(project_root, 'results');
            obj.scripts_dir   = fullfile(project_root, 'scripts_run');
            obj.raw_dir       = fullfile(project_root, 'results_raw');
            obj.cell_def_file = fullfile(obj.inputs_dir, 'cell_def.csv');
            addpath(obj.scripts_dir);
        end

        function run(obj)
            obj.loadAll();
            obj.checkStructure();
            obj.checkVerification();
            obj.printSummary();
        end

        function compute(obj, nworkers, cell_def_file)
            if nargin < 2 || isempty(nworkers),      nworkers = 20; end
            if nargin < 3 || isempty(cell_def_file), cell_def_file = obj.cell_def_file; end
            if nworkers ~= 20
                error('VerifyTriangleInequalities:UnifiedWorkerCount', ...
                    ['The authoritative ver10 wrapper fixes workers=20.  ', ...
                     'Call run_omega_mid_unified_parallel directly for ', ...
                     'a diagnostic run with another worker count.']);
            end
            if ~strcmp(local_canonical_path(cell_def_file), ...
                    local_canonical_path(obj.cell_def_file))
                error('VerifyTriangleInequalities:UnifiedInputFile', ...
                    ['The authoritative publication wrapper uses the ', ...
                     'tracked inputs/cell_def.csv.  Call the unified ', ...
                     'parallel driver directly for another input.']);
            end
            run(fullfile(obj.scripts_dir,'run_omega_mid_unified_full.m'));
        end

        function aggregate(obj)
            if ~exist(obj.results_dir, 'dir'), mkdir(obj.results_dir); end
            for conj_c = {'J1','J2'}
                conj = conj_c{1};
                fname = sprintf('%s_OmegaMid.csv', conj);
                out = fullfile(obj.results_dir, fname);
                rows_map = containers.Map('KeyType','double','ValueType','any');
                workers = dir(fullfile(obj.raw_dir, 'worker_*'));
                for w = 1:numel(workers)
                    p = fullfile(workers(w).folder, workers(w).name, fname);
                    if ~exist(p, 'file'), continue; end
                    fid = fopen(p,'r'); fgetl(fid);
                    while true
                        ln = fgetl(fid); if ~ischar(ln), break; end
                        if isempty(strtrim(ln)), continue; end
                        parts = strsplit(ln, ',');
                        cid = str2double(parts{2});
                        rows_map(cid) = ln;
                    end
                    fclose(fid);
                end
                fid = fopen(out, 'w');
                fprintf(fid, 'conjecture,cell_id,verified,J_lower,status,note,run_timestamp\n');
                cids = sort(cell2mat(rows_map.keys()));
                for k = 1:numel(cids)
                    fprintf(fid, '%s\n', rows_map(cids(k)));
                end
                fclose(fid);
                fprintf('[aggregate] %s: %d rows -> %s\n', conj, numel(cids), out);
            end
        end

        function loadAll(obj)
            J1_path = fullfile(obj.results_dir, 'J1_OmegaMid.csv');
            J2_path = fullfile(obj.results_dir, 'J2_OmegaMid.csv');
            fprintf('[load] %s\n', obj.cell_def_file);
            obj.cell_def = readtable(obj.cell_def_file);
            fprintf('[load] %s\n', J1_path); obj.J1 = readtable(J1_path);
            fprintf('[load] %s\n', J2_path); obj.J2 = readtable(J2_path);
            fprintf('  cell_def: %d cells; J1: %d rows; J2: %d rows\n', ...
                height(obj.cell_def), height(obj.J1), height(obj.J2));
        end

        function checkStructure(obj)
            fprintf('\n[check] structural integrity\n');
            assert(height(obj.cell_def) > 0, 'cell_def empty');
            req_def = {'i','x_inf','x_sup','theta_inf','theta_sup'};
            assert(all(ismember(req_def, obj.cell_def.Properties.VariableNames)), 'cell_def missing columns');
            req_res = {'cell_id','verified','J_lower','status'};
            assert(all(ismember(req_res, obj.J1.Properties.VariableNames)), 'J1 missing columns');
            assert(all(ismember(req_res, obj.J2.Properties.VariableNames)), 'J2 missing columns');
            assert(numel(unique(obj.J1.cell_id))==height(obj.J1), 'J1 cell_id duplicates');
            assert(numel(unique(obj.J2.cell_id))==height(obj.J2), 'J2 cell_id duplicates');
            expected_ids = sort(double(obj.cell_def.i(:)));
            assert(height(obj.J1)==numel(expected_ids), ...
                'J1 row count does not match cell_def');
            assert(height(obj.J2)==numel(expected_ids), ...
                'J2 row count does not match cell_def');
            assert(isequal(sort(double(obj.J1.cell_id(:))),expected_ids), ...
                'J1 cell_id set does not match cell_def.i');
            assert(isequal(sort(double(obj.J2.cell_id(:))),expected_ids), ...
                'J2 cell_id set does not match cell_def.i');
            fprintf(['  OK: columns present, no duplicates, and both ', ...
                'cell-id sets exactly match cell_def.\n']);
        end

        function checkVerification(obj)
            fprintf('\n[check] verified==1 everywhere\n');
            n1u = sum(obj.J1.verified ~= 1); n2u = sum(obj.J2.verified ~= 1);
            if n1u==0, fprintf('  J1: %d/%d rows verified.\n', height(obj.J1), height(obj.J1));
            else, warning('J1 has %d unverified rows', n1u); end
            if n2u==0, fprintf('  J2: %d/%d rows verified.\n', height(obj.J2), height(obj.J2));
            else, warning('J2 has %d unverified rows', n2u); end
            fprintf('\n[check] J_lower finite and > 0 everywhere\n');
            bad1 = ~isfinite(obj.J1.J_lower) | obj.J1.J_lower <= 0;
            bad2 = ~isfinite(obj.J2.J_lower) | obj.J2.J_lower <= 0;
            if ~any(bad1)
                fprintf('  J1: all J_lower finite and positive (min=%.4e).\n', ...
                    min(obj.J1.J_lower));
            else
                warning('J1 has %d nonfinite or nonpositive J_lower rows',sum(bad1));
            end
            if ~any(bad2)
                fprintf('  J2: all J_lower finite and positive (min=%.4e).\n', ...
                    min(obj.J2.J_lower));
            else
                warning('J2 has %d nonfinite or nonpositive J_lower rows',sum(bad2));
            end
            bad_status1 = ~strcmpi(strtrim(string(obj.J1.status)),'ok');
            bad_status2 = ~strcmpi(strtrim(string(obj.J2.status)),'ok');
            if any(bad_status1), warning('J1 has %d non-ok status rows',sum(bad_status1)); end
            if any(bad_status2), warning('J2 has %d non-ok status rows',sum(bad_status2)); end
        end

        function printSummary(obj)
            fprintf('\n========== VERIFICATION SUMMARY ==========\n');
            for conj_c = {'J1','J2'}
                conj = conj_c{1}; T = obj.(conj);
                certified = T.verified == 1 ...
                    & isfinite(T.J_lower) & T.J_lower > 0 ...
                    & strcmpi(strtrim(string(T.status)),'ok');
                n = height(T); v = sum(certified);
                finite_values = T.J_lower(isfinite(T.J_lower));
                if isempty(finite_values)
                    Jmin = NaN; Jmed = NaN; Jmax = NaN;
                else
                    Jmin = min(finite_values);
                    Jmed = median(finite_values);
                    Jmax = max(finite_values);
                end
                fprintf('%s: %d cells  verified=%d (%.3f%%)  J_lower min/med/max = %.3e / %.3e / %.3e\n', ...
                    conj, n, v, v/n*100, Jmin, Jmed, Jmax);
            end
            d1 = sum(contains(string(obj.J1.note), 'derived from subs'));
            d2 = sum(contains(string(obj.J2.note), 'derived from subs'));
            fprintf('Parents verified via subdivision: J1=%d, J2=%d\n', d1, d2);
            certified1 = obj.J1.verified==1 ...
                & isfinite(obj.J1.J_lower) & obj.J1.J_lower > 0 ...
                & strcmpi(strtrim(string(obj.J1.status)),'ok');
            certified2 = obj.J2.verified==1 ...
                & isfinite(obj.J2.J_lower) & obj.J2.J_lower > 0 ...
                & strcmpi(strtrim(string(obj.J2.status)),'ok');
            all_ok = all(certified1) && all(certified2);
            if all_ok
                fprintf('\n*** VERIFIED: J1 >= 0 and J2 >= 0 on all cells of Omega_mid. ***\n');
            else
                fprintf('\n!!! INCOMPLETE: some cells unverified or J_lower <= 0.\n');
                error('VerifyTriangleInequalities:IncompleteCertificate', ...
                    ['Committed Omega_mid results are incomplete: every ', ...
                     'row must have verified=1, finite J_lower>0, and status=ok.']);
            end
        end
    end
end


function path_value = local_canonical_path(path_value)
path_value = char(java.io.File(char(path_value)).getCanonicalPath());
end
