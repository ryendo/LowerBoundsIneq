function my_intlab_config()
%MY_INTLAB_CONFIG  Initialize INTLAB and add all project paths.
%   Set INTLAB_ROOT, GMSH_COMMAND, and LOWERBOUNDS_MESH_PATH in the
%   environment when needed, then call this function once at the top of
%   any driver script.  The constants below are portable fallbacks.

    % ============================================================
    % User-configurable constants
    % ============================================================
    % Path to the gmsh executable (tested with gmsh 4.8.4).
    % If your gmsh is on PATH as a regular command, 'gmsh' also works.
    GMSH_COMMAND = '/usr/bin/gmsh';

    % (optional) Folder for temporary mesh .geo/.msh files.
    % Empty => use $LOWERBOUNDS_MESH_PATH when set, otherwise a private
    % directory below tempdir.  Keeping generated meshes out of src/ avoids
    % dirtying a checkout and prevents parallel workers from sharing files.
    MESH_PATH    = '';
    % ============================================================

    project_root = fileparts(mfilename('fullpath'));

    global INTERVAL_MODE gmsh_command mesh_path
    % Fail closed: interval mode is enabled only after STARTINTLAB returns
    % successfully in this process.
    INTERVAL_MODE = 0;
    gmsh_command = getenv('GMSH_COMMAND');
    if isempty(gmsh_command)
        gmsh_command = GMSH_COMMAND;
    end
    if isempty(MESH_PATH)
        mesh_path = getenv('LOWERBOUNDS_MESH_PATH');
    else
        mesh_path = MESH_PATH;
    end
    if isempty(mesh_path)
        mesh_path = fullfile(tempdir, 'LowerBoundsIneq-mesh');
    end
    if ~exist(mesh_path, 'dir')
        mkdir(mesh_path);
    end

    % INTLAB is intentionally not vendored.  Prefer an explicit environment
    % variable on shared/HPC installations and retain the historical local
    % folder as a fallback for existing users.
    intlab_root = getenv('INTLAB_ROOT');
    if isempty(intlab_root)
        intlab_root = fullfile(project_root, 'Intlab_V12');
    end
    if ~isfolder(intlab_root)
        error(['INTLAB directory not found. Set INTLAB_ROOT or place ', ...
               'Intlab_V12 below the project root.']);
    end

    % Add project source directories
    addpath(intlab_root);
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

    % STARTINTLAB loads and then rewrites a shared MAT file containing
    % proof-critical INTLAB_CONST data.  Serialize this section across
    % clients and parpool workers so no process can observe a partial save.
    intlab_lock = ver10_acquire_intlab_init_lock(intlab_root);
    try
        evalc('startintlab');
        INTERVAL_MODE = 1;
        clear intlab_lock
        fprintf('INTLAB initialized.\n');
        fprintf('Project root: %s\n', project_root);
        fprintf('INTLAB root: %s\n', intlab_root);
        fprintf('Mesh workspace: %s\n', mesh_path);
    catch ME
        error('startintlab failed: %s', ME.message);
    end
end
