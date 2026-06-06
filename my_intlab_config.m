function my_intlab_config()
%MY_INTLAB_CONFIG  Initialize INTLAB and add all project paths.
%   Edit the three constants below to match your local system, then call
%   this function once at the top of any driver script.

    % ============================================================
    % User-configurable constants
    % ============================================================
    % Path to the gmsh executable (tested with gmsh 4.8.4).
    % If your gmsh is on PATH as a regular command, 'gmsh' also works.
    GMSH_COMMAND = '/usr/bin/gmsh';

    % (optional) Folder for temporary mesh .geo/.msh files.
    % Empty => use the project's src/mesh/ directory.
    MESH_PATH    = '';
    % ============================================================

    project_root = fileparts(mfilename('fullpath'));

    global INTERVAL_MODE gmsh_command mesh_path
    INTERVAL_MODE = 1;
    gmsh_command  = GMSH_COMMAND;
    if isempty(MESH_PATH)
        mesh_path = fullfile(project_root, 'src', 'mesh');
    else
        mesh_path = MESH_PATH;
    end

    % Add project source directories
    addpath(fullfile(project_root, 'Intlab_V12'));
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

    try
        evalc('startintlab');
        fprintf('INTLAB initialized.\n');
        fprintf('Project root: %s\n', project_root);
    catch ME
        error('startintlab failed: %s', ME.message);
    end
end
