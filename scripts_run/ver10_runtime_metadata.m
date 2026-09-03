function runtime = ver10_runtime_metadata()
%VER10_RUNTIME_METADATA Record portable dependency/runtime identifiers.

runtime = struct();
if exist('OCTAVE_VERSION','builtin')
    runtime.engine = 'GNU Octave';
    runtime.engine_version = OCTAVE_VERSION;
    runtime.release = '';
else
    runtime.engine = 'MATLAB';
    runtime.engine_version = version;
    runtime.release = version('-release');
end
runtime.computer = computer;
runtime.interval_toolbox = 'INTLAB';
intlab_root = getenv('INTLAB_ROOT');
if isempty(intlab_root)
    runtime.intlab_root_name = 'project-local Intlab_V12 fallback';
else
    [~,runtime.intlab_root_name] = fileparts(intlab_root);
end
startintlab_path = which('startintlab');
intval_constructor_path = which('intval');
runtime.startintlab_file = local_basename(startintlab_path);
runtime.intval_constructor_file = ...
    local_basename(intval_constructor_path);
runtime.startintlab_sha256 = local_file_sha256( ...
    startintlab_path);
runtime.intval_constructor_sha256 = local_file_sha256( ...
    intval_constructor_path);
if isempty(intlab_root) && ~isempty(startintlab_path)
    intlab_root = fileparts(startintlab_path);
end
[runtime.intlab_tree_sha256,runtime.intlab_tree_file_count, ...
    runtime.intlab_tree_total_bytes] = local_directory_sha256(intlab_root);
global gmsh_command
gmsh = gmsh_command;
if isempty(gmsh)
    gmsh = getenv('GMSH_COMMAND');
end
if isempty(gmsh)
    gmsh = 'gmsh';
end
[~,gmsh_name,gmsh_extension] = fileparts(gmsh);
runtime.gmsh_executable = [gmsh_name,gmsh_extension];
[resolve_status,resolved_gmsh] = system( ...
    ['command -v ',local_shell_quote(gmsh)]);
resolved_gmsh = strtrim(resolved_gmsh);
if resolve_status ~= 0 || isempty(resolved_gmsh)
    resolved_gmsh = gmsh;
end
[gmsh_status,gmsh_output] = system( ...
    ['env -u LD_LIBRARY_PATH -u LD_PRELOAD ', ...
     local_shell_quote(resolved_gmsh),' --version']);
runtime.gmsh_version_status = gmsh_status;
if gmsh_status == 0
    runtime.gmsh_version = strtrim(gmsh_output);
else
    runtime.gmsh_version = '';
    runtime.gmsh_version_error = strtrim(gmsh_output);
end
if isfile(resolved_gmsh)
    runtime.gmsh_binary_sha256 = ver10_file_sha256(resolved_gmsh);
else
    runtime.gmsh_binary_sha256 = '';
end
end


function [digest,file_count,total_bytes] = local_directory_sha256(directory)
% Bind the complete trusted INTLAB installation, not only two entry files.
digest = '';
file_count = 0;
total_bytes = 0;
if isempty(directory) || ~isfolder(directory)
    return
end

entries = dir(fullfile(directory,'**','*'));
entries = entries(~[entries.isdir]);
if isempty(entries)
    return
end
full_names = arrayfun(@(entry) fullfile(entry.folder,entry.name), ...
    entries,'UniformOutput',false);
root_prefix = [char(directory),filesep];
relative_names = full_names;
for k = 1:numel(relative_names)
    local_name = relative_names{k};
    if ~startsWith(local_name,root_prefix)
        error('ver10_runtime_metadata:BadIntlabTreeEntry', ...
            'INTLAB tree entry is outside the recorded root: %s',local_name);
    end
    relative_names{k} = strrep( ...
        local_name(numel(root_prefix)+1:end),filesep,'/');
end
[relative_names,order] = sort(relative_names);
full_names = full_names(order);
entries = entries(order);

try
    engine = javaMethod( ...
        'getInstance','java.security.MessageDigest','SHA-256');
catch ME
    error('ver10_runtime_metadata:NoSHA256', ...
        'Cannot initialize the SHA-256 engine: %s',ME.message);
end
local_digest_update(engine,uint8('ver10-directory-sha256-v1'));
local_digest_update(engine,uint8(0));
for k = 1:numel(full_names)
    path_bytes = unicode2native(relative_names{k},'UTF-8');
    size_bytes = uint8(sprintf('%d',entries(k).bytes));
    local_digest_update(engine,path_bytes);
    local_digest_update(engine,uint8(0));
    local_digest_update(engine,size_bytes);
    local_digest_update(engine,uint8(0));

    fid = fopen(full_names{k},'rb');
    if fid < 0
        error('ver10_runtime_metadata:CannotHashIntlabFile', ...
            'Cannot open INTLAB tree file %s.',full_names{k});
    end
    cleanup = onCleanup(@() fclose(fid));
    while true
        bytes = fread(fid,1024*1024,'*uint8');
        if isempty(bytes)
            break
        end
        local_digest_update(engine,bytes);
    end
    clear cleanup
    local_digest_update(engine,uint8(0));
end
raw = typecast(engine.digest(),'uint8');
digest = lower(reshape(dec2hex(raw,2).',1,[]));
file_count = numel(entries);
total_bytes = sum(double([entries.bytes]));
end


function local_digest_update(engine,bytes)
engine.update(typecast(uint8(bytes(:)),'int8'));
end


function digest = local_file_sha256(filename)
if isempty(filename) || ~isfile(filename)
    digest = '';
else
    digest = ver10_file_sha256(filename);
end
end


function name = local_basename(filename)
if isempty(filename)
    name = '';
else
    [~,stem,extension] = fileparts(filename);
    name = [stem,extension];
end
end


function quoted = local_shell_quote(value)
% POSIX single-quote escaping; all production runs use Linux/macOS shells.
quoted = ['''',strrep(char(value),'''','''"''"'''),''''];
end
