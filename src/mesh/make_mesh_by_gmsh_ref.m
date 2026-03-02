function mesh = make_mesh_by_gmsh_ref(a,b,h)
global mesh_path gmsh_command

this_dir = fileparts(mfilename('fullpath'));              % .../src/mesh
if isempty(mesh_path), mesh_path = this_dir; end
sh = fullfile(this_dir,'create_mesh.sh');

V  = I_intval([0 0; 1 0; a b]);
FV = [0 0; 1 0; I_mid(a) I_mid(b)];

fid = fopen(fullfile(mesh_path,'temp.geo'),'w');
for i=1:3
    fprintf(fid,'Point(%d)={%.17g,%.17g,0,%.17g};\n',i,FV(i,1),FV(i,2),I_mid(h));
end
fprintf(fid,'Line(1)={1,2};Line(2)={2,3};Line(3)={3,1};\n');
fprintf(fid,'Line Loop(1)={1,2,3};Plane Surface(1)={1};\n');
fclose(fid);

cmd = sprintf('bash "%s" "%s" "%s"', sh, gmsh_command, mesh_path);
[st,out] = system(cmd); if st~=0, error("Mesh generation failed:\n%s",out); end


mesh = gmshread(fullfile(mesh_path,'temp.msh'));
mesh.domain = V;
mesh = apply_exact_boundary_point_setting(mesh,V);
end
