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
intlab_data_relative_names = local_intlab_data_relative_names(intlab_root);
[runtime.intlab_tree_sha256,runtime.intlab_tree_file_count, ...
    runtime.intlab_tree_total_bytes,normalized_records] = ...
    ver10_directory_sha256(intlab_root,intlab_data_relative_names, ...
        @() ver10_acquire_intlab_init_lock(intlab_root));
runtime.intlab_tree_hash_policy = [ ...
    'recursive-path-size-content-sha256-v2-with-exact-intlab-data-', ...
    'mat-v5-created-on-timestamp-normalization'];
runtime.intlab_data_mat_records = normalized_records;
runtime.intlab_data_mat_write_serialization = ...
    'java-nio-advisory-lock-adjacent-to-intlab-root-v1';
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


function relative_names = local_intlab_data_relative_names(intlab_root)
% Identify the exact platform-data MAT written and loaded by INTLAB.  Its
% serialized INTLAB_CONST payload is proof-critical and remains hashed;
% only the timestamp-bearing MAT-v5 description is normalized downstream.
relative_names = {};
if isempty(intlab_root) || ~isfolder(intlab_root)
    return
end

candidate = '';
if exist('intvalinit','file') == 2
    try
        candidate = char(intvalinit('intlabdata'));
    catch ME
        error('ver10_runtime_metadata:IntlabDataLookupFailed', ...
            'INTLAB data MAT lookup failed: %s',ME.message);
    end
end
if isempty(candidate)
    if exist('OCTAVE_VERSION','builtin')
        prefix = ['Octave_',version,'_Intlab_Version_'];
    else
        prefix = ['Matlab_',version,'_Intlab_Version_'];
    end
    entries = dir(fullfile(intlab_root,[prefix,'*.mat']));
    entries = entries(~[entries.isdir]);
    if numel(entries) == 1
        candidate = fullfile(entries(1).folder,entries(1).name);
    elseif numel(entries) > 1
        error('ver10_runtime_metadata:AmbiguousIntlabData', ...
            'Multiple current-runtime INTLAB data MAT files were found.');
    end
end
if isempty(candidate)
    return
end

[candidate_parent,stem,extension] = fileparts(candidate);
candidate_name = [stem,extension];
if exist('OCTAVE_VERSION','builtin')
    exact_prefix = ['Octave_',version,'_Intlab_Version_'];
else
    exact_prefix = ['Matlab_',version,'_Intlab_Version_'];
end
matching = dir(fullfile(intlab_root,[exact_prefix,'*.mat']));
matching = matching(~[matching.isdir]);
if ~local_same_directory(candidate_parent,intlab_root) ...
        || ~local_is_unsymlinked_root_child(candidate,intlab_root, ...
            candidate_name) ...
        || ~startsWith(candidate_name,exact_prefix) ...
        || ~endsWith(candidate_name,'.mat') ...
        || numel(matching) ~= 1 ...
        || ~strcmp(matching(1).name,candidate_name) ...
        || ~isfile(candidate)
    error('ver10_runtime_metadata:BadIntlabDataPath', ...
        ['INTLAB reported an invalid platform-data MAT path; refusing ', ...
         'to normalize any tree entry.']);
end
relative_names = {candidate_name};
end


function tf = local_same_directory(first,second)
try
    first_file = javaObject('java.io.File',char(first));
    second_file = javaObject('java.io.File',char(second));
    tf = strcmp(char(first_file.getCanonicalPath()), ...
        char(second_file.getCanonicalPath()));
catch
    tf = strcmp(char(first),char(second));
end
end


function tf = local_is_unsymlinked_root_child(candidate,root,name)
try
    root_file = javaObject('java.io.File',char(root));
    candidate_file = javaObject('java.io.File',char(candidate));
    expected = fullfile(char(root_file.getCanonicalPath()),char(name));
    tf = candidate_file.isFile() ...
        && strcmp(char(candidate_file.getCanonicalPath()),expected);
catch
    tf = false;
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
