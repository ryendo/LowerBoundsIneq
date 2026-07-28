function omega_up_all_prepare_worker(project_root,mode)
%OMEGA_UP_ALL_PREPARE_WORKER  Initialize one serial/parpool runtime.

mode = lower(char(mode));
addpath(project_root);
addpath(fullfile(project_root,'scripts_run'));

if strcmp(mode,'interval')
    my_intlab_config_worker();
    return;
end
if ~strcmp(mode,'double')
    error('omega_up_all_prepare_worker:BadMode', ...
        'mode must be ''interval'' or ''double''.');
end

addpath(fullfile(project_root,'src','algorithms'));
addpath(fullfile(project_root,'src','fem'));
addpath(fullfile(project_root,'src','mesh'));
addpath(fullfile(project_root,'src','interval'));
addpath(fullfile(project_root,'src'));
addpath(fullfile(project_root,'src','lib','VFEM2D', ...
    'lib_eigenvalue_bound'));
addpath(fullfile(project_root,'src','lib','VFEM2D_revised'));
addpath(fullfile(project_root,'src','lib','veigs'));
addpath(fullfile(project_root,'tests'));

global INTERVAL_MODE gmsh_command mesh_path
INTERVAL_MODE = 0;
gmsh_from_environment = getenv('GMSH_COMMAND');
if isempty(gmsh_from_environment)
    gmsh_from_environment = '/usr/bin/gmsh';
end
gmsh_command = gmsh_from_environment;
mesh_path = fullfile(tempdir,sprintf( ...
    'LowerBoundsIneq-mesh-pid-%d',local_process_id()));
if ~exist(mesh_path,'dir')
    mkdir(mesh_path);
end
end


function pid = local_process_id()
pid = 0;
try
    pid = feature('getpid');
catch
    try
        pid = getpid();
    catch
        pid = round(1e9*rand());
    end
end
end
