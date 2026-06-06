%% Main Script: Visualize Verification Results
clear; clc; close all;

try
    Cell_Def = readtable('cell_def.csv');
    if ismember('i', Cell_Def.Properties.VariableNames)
        Cell_Def = renamevars(Cell_Def, 'i', 'cell_id');
    end
catch ME
    error('Error loading CSV files. Ensure files exist in the current folder.\n%s', ME.message);
end

visualize_cells(Cell_Def)

%% --- Local Function: Visualization Logic ---
function visualize_cells(Cell_Def)
    Data = Cell_Def;
    num_cells = height(Data);
    
    x_inf = Data.x_inf;
    x_sup = Data.x_sup;
    t_inf = Data.theta_inf;
    t_sup = Data.theta_sup;
    
    face_colors = repmat([0.6 0.8 1.0], num_cells, 1);
    
    X_log = [x_inf'; x_sup'; x_sup'; x_inf'];
    T_log = [t_inf'; t_inf'; t_sup'; t_sup'];
    
    X_phys = X_log;
    Y_phys = X_log .* tan(T_log);
    
    dark_red = [0.6 0 0];
    
    % --- Plot 1: (x, theta) Domain ---
    figure('Name', 'Cells: (x, theta)', 'Color', 'w');
    p1 = patch(X_log, T_log, 'w');
    
    if exist('set_patch_properties', 'file')
        set_patch_properties(p1, face_colors);
    else
        set(p1, 'FaceColor', 'flat', 'FaceVertexCData', face_colors, 'EdgeColor', 'k');
    end
    
    hold on;
    % 領域の x の範囲は 0.5 から 1 になる
    x_bnd = linspace(0.5, 1, 100);
    
    theta_arc = acos(x_bnd); 
    theta_line = zeros(size(x_bnd)); % y=0 のため theta=0
    
    plot(x_bnd, theta_arc, 'Color', dark_red, 'LineWidth', 2); 
    plot(x_bnd, theta_line, 'Color', dark_red, 'LineWidth', 2); 
    plot([0.5, 0.5], [0, acos(0.5)], 'Color', dark_red, 'LineWidth', 2); 
    
    xlabel('x');
    ylabel('\theta');
    axis tight; 
    grid on;
    hold off;
    
    % --- Plot 2: (x, y) Domain ---
    figure('Name', 'Cells: (x, y)', 'Color', 'w');
    p2 = patch(X_phys, Y_phys, 'w');
    
    if exist('set_patch_properties', 'file')
        set_patch_properties(p2, face_colors);
    else
        set(p2, 'FaceColor', 'flat', 'FaceVertexCData', face_colors, 'EdgeColor', 'k');
    end
    
    hold on;
    y_arc = sqrt(1 - x_bnd.^2);
    y_line = zeros(size(x_bnd));
    
    plot(x_bnd, y_arc, 'Color', dark_red, 'LineWidth', 2); 
    plot(x_bnd, y_line, 'Color', dark_red, 'LineWidth', 2); 
    plot([0.5, 0.5], [0, sqrt(1 - 0.5^2)], 'Color', dark_red, 'LineWidth', 2); 
    
    xlabel('x');
    ylabel('y');
    axis equal; 
    axis tight;
    grid on;
    hold off;
end