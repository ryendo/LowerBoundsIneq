classdef VerifyTriangleInequalities < handle
%VERIFYTRIANGLEINEQUALITIES  Main class for computer-assisted verification of
%   sharp Dirichlet Laplacian eigenvalue inequalities on planar triangles.
%
% Usage:
%   (1) Validate committed results without recomputing:
%           v = VerifyTriangleInequalities();
%           v.run();
%   (2) Reproduce the full computation from scratch:
%           v = VerifyTriangleInequalities();
%           v.compute(20);  % use 20 parallel workers
%
% .compute() runs parallel FEM verification over inputs/cell_def.csv,
% AUTO-AGGREGATES per-worker outputs into results/J1_OmegaMid.csv and
% results/J2_OmegaMid.csv, and calls .run() to validate everything passes.
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
            cd(obj.project_root);
            my_intlab_config();
            t0 = tic;
            fprintf('\n=== VerifyTriangleInequalities.compute ===\n');
            fprintf('cell_def : %s\n', cell_def_file);
            fprintf('workers  : %d\n', nworkers);
            if ~exist(obj.raw_dir,'dir'), mkdir(obj.raw_dir); end
            fprintf('\n--- J1 ---\n');
            run_parallel_omegamid('J1', [], nworkers, cell_def_file, obj.raw_dir);
            fprintf('\n--- J2 ---\n');
            run_parallel_omegamid('J2', [], nworkers, cell_def_file, obj.raw_dir);
            fprintf('\n--- aggregate raw -> results/ ---\n');
            obj.aggregate();
            fprintf('\n--- validate ---\n');
            obj.run();
            fprintf('\n=== total wall time: %.2f hours ===\n', toc(t0)/3600);
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
            fprintf('  OK: columns present, no duplicates.\n');
        end

        function checkVerification(obj)
            fprintf('\n[check] verified==1 everywhere\n');
            n1u = sum(obj.J1.verified ~= 1); n2u = sum(obj.J2.verified ~= 1);
            if n1u==0, fprintf('  J1: %d/%d rows verified.\n', height(obj.J1), height(obj.J1));
            else, warning('J1 has %d unverified rows', n1u); end
            if n2u==0, fprintf('  J2: %d/%d rows verified.\n', height(obj.J2), height(obj.J2));
            else, warning('J2 has %d unverified rows', n2u); end
            fprintf('\n[check] J_lower > 0 everywhere\n');
            n1n = sum(obj.J1.J_lower <= 0); n2n = sum(obj.J2.J_lower <= 0);
            if n1n==0, fprintf('  J1: all J_lower positive (min=%.4e).\n', min(obj.J1.J_lower));
            else, warning('J1 has %d rows with J_lower <= 0', n1n); end
            if n2n==0, fprintf('  J2: all J_lower positive (min=%.4e).\n', min(obj.J2.J_lower));
            else, warning('J2 has %d rows with J_lower <= 0', n2n); end
        end

        function printSummary(obj)
            fprintf('\n========== VERIFICATION SUMMARY ==========\n');
            for conj_c = {'J1','J2'}
                conj = conj_c{1}; T = obj.(conj);
                n = height(T); v = sum(T.verified == 1);
                Jmin = min(T.J_lower); Jmed = median(T.J_lower); Jmax = max(T.J_lower);
                fprintf('%s: %d cells  verified=%d (%.3f%%)  J_lower min/med/max = %.3e / %.3e / %.3e\n', ...
                    conj, n, v, v/n*100, Jmin, Jmed, Jmax);
            end
            d1 = sum(contains(string(obj.J1.note), 'derived from subs'));
            d2 = sum(contains(string(obj.J2.note), 'derived from subs'));
            fprintf('Parents verified via subdivision: J1=%d, J2=%d\n', d1, d2);
            all_ok = (sum(obj.J1.verified==1)==height(obj.J1)) && (sum(obj.J2.verified==1)==height(obj.J2)) && ...
                     all(obj.J1.J_lower > 0) && all(obj.J2.J_lower > 0);
            if all_ok
                fprintf('\n*** VERIFIED: J1 >= 0 and J2 >= 0 on all cells of Omega_mid. ***\n');
            else
                fprintf('\n!!! INCOMPLETE: some cells unverified or J_lower <= 0.\n');
            end
        end
    end
end
