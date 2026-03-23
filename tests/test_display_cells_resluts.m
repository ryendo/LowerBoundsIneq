%% Main Script: Visualize Verification Results for J1 and J2
clear; clc; close all;

cell_def_path = 'inputs/cell_def.csv';
results_path = 'results/J1_OmegaMid.csv';

% --- 1. Load Data ---
try
    % Read the geometry definition
    optsCell = detectImportOptions(cell_def_path, 'Delimiter', ',');
    Cell_Def = readtable(cell_def_path, optsCell);
    
    % Read the conjecture results
    optsJ1 = detectImportOptions(results_path, 'Delimiter', ',');
    optsJ1 = setvartype(optsJ1, {'conjecture','status','note','run_timestamp'}, 'string');
    J1_Data = readtable(results_path, optsJ1);
    
    % Rename 'i' to 'cell_id' in Cell_Def to match J tables for joining
    if ismember('i', Cell_Def.Properties.VariableNames)
        Cell_Def = renamevars(Cell_Def, 'i', 'cell_id');
    end

    % Make sure join keys have the same type
    J1_Data.cell_id = double(J1_Data.cell_id);
    Cell_Def.cell_id = double(Cell_Def.cell_id);
    
catch ME
    error('Error loading CSV files. Ensure files exist in the current folder.\n%s', ME.message);
end

%% --- 2. Visualize J1 ---
visualize_conjecture(J1_Data, Cell_Def, 'J1');

%% --- 3. Visualize J2 ---
% visualize_conjecture(J2_Data, Cell_Def, 'J2');


%% --- Local Function: Visualization Logic ---
function visualize_conjecture(J_Data, Cell_Def, label_name)
    % Join the verification data with the cell definitions
    % 'innerjoin' ensures we only plot cells present in both files

    Data = innerjoin(J_Data, Cell_Def, 'Keys', 'cell_id');
    
    num_cells = height(Data);
    
    % Extract boundaries
    x_inf = Data.x_inf;
    x_sup = Data.x_sup;
    t_inf = Data.theta_inf;
    t_sup = Data.theta_sup;
    verified = Data.verified;
    
    % --- Prepare Colors ---
    % Green [0 1 0] for verified (1)
    % Red   [1 0 0] for unverified (0)
    face_colors = zeros(num_cells, 3);
    face_colors(verified == 1, :) = repmat([0 1 0], sum(verified == 1), 1);
    face_colors(verified == 0, :) = repmat([1, 0, 0], sum(verified == 0), 1);
    
    % --- Prepare Vertices for Patch (Counter-Clockwise) ---
    % Rows: 4 vertices, Cols: N cells
    
    % 1. Logical Domain (x, theta)
    X_log = [x_inf'; x_sup'; x_sup'; x_inf'];
    T_log = [t_inf'; t_inf'; t_sup'; t_sup'];
    
    % 2. Physical Domain (x, y)
    % Transformation: y = x * tan(theta)
    X_phys = X_log;
    Y_phys = X_log .* tan(T_log);
    
    % --- Plot 1: (x, theta) Domain ---
    figure('Name', [label_name ': (x, theta)'], 'Color', 'w');
    p1 = patch(X_log, T_log, 'w');
    set_patch_properties(p1, face_colors);
    
    title([label_name ': Verification on (x, \theta) plane']);
    xlabel('x');
    ylabel('\theta');
    axis tight; 
    grid on;
    add_custom_legend(verified);
    
    % --- Plot 2: (x, y) Domain ---
    figure('Name', [label_name ': (x, y)'], 'Color', 'w');
    p2 = patch(X_phys, Y_phys, 'w');
    set_patch_properties(p2, face_colors);
    
    title([label_name ': Verification on (x, y) plane']);
    xlabel('x');
    ylabel('y');
    axis equal; % Important to maintain physical aspect ratio
    axis tight;
    grid on;
    add_custom_legend(verified);
end

function set_patch_properties(p, colors)
    % Helper to apply colors and style to the patch object
    p.FaceVertexCData = colors;
    p.FaceColor = 'flat';
    p.EdgeColor = 'none'; % Turn off edges to prevent black clutter
    % p.EdgeColor = 'k'; p.EdgeAlpha = 0.1; % Uncomment for faint grid lines
end

function add_custom_legend(verified_array)
    % Adds a manual legend for Green/Red status
    hold on;
    h_ver = plot(nan, nan, 's', 'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', 'none');
    h_unver = plot(nan, nan, 's', 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none');
    
    legend_entries = [];
    legend_labels = {};
    
    if any(verified_array == 1)
        legend_entries = [legend_entries, h_ver];
        legend_labels{end+1} = 'Verified (1)';
    end
    if any(verified_array == 0)
        legend_entries = [legend_entries, h_unver];
        legend_labels{end+1} = 'Unverified (0)';
    end
    
    if ~isempty(legend_entries)
        legend(legend_entries, legend_labels, 'Location', 'best');
    end
    hold off;
end