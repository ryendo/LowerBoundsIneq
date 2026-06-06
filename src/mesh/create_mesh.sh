#!/usr/bin/bash
# args: gmsh_command mesh_path suffix
# Generates a 2D mesh from <mesh_path>/temp<suffix>.geo.

gmsh_command=$1
mesh_path=$2
suffix=$3

# Force gmsh to use the system C++ runtime rather than the one MATLAB
# injects via LD_LIBRARY_PATH; this is needed on Linux hosts where the
# MATLAB-provided libstdc++ is older than the one gmsh was built against.
export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6
unset LD_LIBRARY_PATH

# If the gmsh binary is a Python-wrapped script (some conda distributions),
# ensure a python interpreter is on PATH. Uncomment and edit if needed:
# export PATH=$HOME/miniconda3/bin:$PATH

${gmsh_command} ${mesh_path}/temp${suffix}.geo -2 -format msh2 -o ${mesh_path}/temp${suffix}.msh \
    -setnumber Mesh.Binary 0 \
    -setnumber Mesh.ElementOrder 1
