%RUN_MUHAT1_RAYLEIGH_HPC  HPC/INTLAB certificate for t=0.38.

project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(project_root, 'src', 'degenerate'));

result = muhat1_rayleigh_cover( ...
    'mode', 'interval', ...
    't', '0.38', ...
    'degree', 10, ...
    's_cells', 200, ...
    'down_y_max', 0.19, ...
    'down_y_cells', 40, ...
    'progress_every', 10);

output_dir = fullfile(project_root, 'results', ...
    'muhat1_rayleigh_t038_hpc_20260730');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

cell_table = struct2table(result.cells);
writetable(cell_table, fullfile(output_dir, 'cells.csv'));

summary = rmfield(result, 'cells');
summary.matlab_version = version;
summary.computer = computer;
summary.intlab_root = getenv('INTLAB_ROOT');
summary.generated_utc = char(datetime('now', ...
    'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z'''));

fid = fopen(fullfile(output_dir, 'summary.json'), 'w');
if fid < 0
    error('Could not open summary.json for writing.');
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', jsonencode(summary, 'PrettyPrint', true));

fprintf('t                         : %.17g\n', summary.t);
fprintf('degree                    : %d\n', summary.degree);
fprintf('s cells                   : %d\n', summary.s_cells);
fprintf('max certified Rayleigh UB : %.17g\n', summary.max_rayleigh_upper);
fprintf('max cell                  : [%g,%g]\n', ...
    summary.max_cell(1), summary.max_cell(2));
fprintf('target 11.5 verified      : %d\n', summary.target_11_5_verified);
fprintf('max down scaled upper     : %.17g\n', ...
    summary.max_down_scaled_upper);
fprintf('down conjecture verified  : %d\n', ...
    summary.down_conjecture_verified);
