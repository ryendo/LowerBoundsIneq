function portable = ver10_portable_path(project_root,filename)
%VER10_PORTABLE_PATH Store project-relative result paths in manifests.

project_root = char(project_root);
filename = char(filename);
prefix = [project_root,filesep];
if startsWith(filename,prefix)
    portable = filename((length(prefix)+1):end);
else
    [~,name,extension] = fileparts(filename);
    portable = fullfile('external-output',[name,extension]);
end
portable = strrep(portable,filesep,'/');
end
