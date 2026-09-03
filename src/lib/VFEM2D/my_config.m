%The path of the codes for switch between verified computing and approximate computing.
addpath('./lib_eigenvalue_bound/');
addpath('./lib_mesh/');
addpath('./mode_switch_interface/');  
addpath("./veigs/");

%The path of INTLAB toolbox and initialization.
addpath('~/app/Intlab_V12/');

clear INTERVAL_MODE;
global INTERVAL_MODE;

%INTERVAL_MODE=1; for rigorous computing based on interval arithmetic.
%INTERVAL_MODE=0; for approximate computing with rounding error inside.
INTERVAL_MODE=1;

if INTERVAL_MODE == 1
    project_root = fileparts(fileparts(fileparts( ...
        fileparts(mfilename('fullpath')))));
    addpath(fullfile(project_root,'scripts_run'));
    intlab_root = fileparts(which('startintlab'));
    intlab_lock = ver10_acquire_intlab_init_lock(intlab_root);
    try
        startintlab;
        clear intlab_lock
    catch ME
        clear intlab_lock
        rethrow(ME);
    end
end
