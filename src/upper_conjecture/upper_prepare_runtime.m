function upper_prepare_runtime(project_root, mode)
%UPPER_PREPARE_RUNTIME Add verifier paths and initialize INTLAB on a worker.
persistent interval_ready

addpath(fullfile(project_root, 'src', 'upper_conjecture'));
addpath(fullfile(project_root, 'src', 'interval'));
addpath(fullfile(project_root, 'src', 'mesh'));
addpath(fullfile(project_root, 'src', 'lib', 'VFEM2D', 'lib_eigenvalue_bound'));
addpath(fullfile(project_root, 'src', 'lib', 'VFEM2D_revised'));
addpath(fullfile(project_root, 'src', 'lib', 'veigs'));
addpath(fullfile(project_root, 'scripts_run'));

global INTERVAL_MODE gmsh_command mesh_path
INTERVAL_MODE = 0;
if isempty(gmsh_command)
    gmsh_command = getenv('GMSH_COMMAND');
    if isempty(gmsh_command)
        gmsh_command = 'gmsh';
    end
end
if isempty(mesh_path)
    mesh_path = fullfile(tempdir,'LowerBoundsIneq-upper-mesh');
end
if ~isfolder(mesh_path)
    mkdir(mesh_path);
end
if ~strcmpi(mode, 'interval')
    return;
end
if ~isempty(interval_ready) && interval_ready
    INTERVAL_MODE = 1;
    return;
end

intlab_root = getenv('INTLAB_ROOT');
if isempty(intlab_root)
    intlab_root = fullfile(project_root, 'Intlab_V12');
end
if ~isfolder(intlab_root)
    error(['INTLAB not found. Set INTLAB_ROOT to an INTLAB installation ', ...
           '(for example /path/to/Intlab_V12).']);
end
addpath(intlab_root);
intlab_lock = ver10_acquire_intlab_init_lock(intlab_root);
evalc('startintlab');
clear intlab_lock
interval_ready = true;
INTERVAL_MODE = 1;
end
