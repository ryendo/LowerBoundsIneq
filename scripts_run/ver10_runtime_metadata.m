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
