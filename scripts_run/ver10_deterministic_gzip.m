function certificate = ver10_deterministic_gzip(source,destination)
%VER10_DETERMINISTIC_GZIP Create and verify a timestamp-free gzip artifact.
%
% GNU-compatible `gzip -n -9` omits the source name and timestamp.  Two
% independent encodings must be byte-identical, and expanding the result
% must reproduce the source SHA-256 before the destination is published.

source = char(source);
destination = char(destination);
if exist(source,'file') ~= 2
    error('ver10_deterministic_gzip:MissingSource', ...
        'Cannot compress missing source %s.',source);
end
directory = fileparts(destination);
if isempty(directory)
    directory = pwd;
elseif ~exist(directory,'dir')
    mkdir(directory);
end
gzip_command = getenv('GZIP_COMMAND');
if isempty(gzip_command)
    gzip_command = 'gzip';
end

first = [tempname(directory),'.gz'];
second = [tempname(directory),'.gz'];
expanded = [tempname(directory),'.csv'];
cleanup = onCleanup(@() local_cleanup({first,second,expanded}));
local_run_gzip(gzip_command,source,first);
local_run_gzip(gzip_command,source,second);
first_sha = ver10_file_sha256(first);
second_sha = ver10_file_sha256(second);
if ~strcmp(first_sha,second_sha)
    error('ver10_deterministic_gzip:NondeterministicOutput', ...
        'Two gzip -n encodings of %s had different SHA-256 hashes.', ...
        source);
end
local_expand_gzip(gzip_command,first,expanded);
content_sha = ver10_file_sha256(source);
expanded_sha = ver10_file_sha256(expanded);
if ~strcmp(content_sha,expanded_sha)
    error('ver10_deterministic_gzip:ExpansionMismatch', ...
        'Expanding the gzip artifact did not reproduce %s.',source);
end
[ok,message] = movefile(first,destination,'f');
if ~ok
    error('ver10_deterministic_gzip:PublishFailed', ...
        'Could not publish %s: %s',destination,message);
end

certificate = struct( ...
    'schema','lowerboundsineq.deterministic_gzip.v1', ...
    'command','gzip -n -9', ...
    'file_sha256',first_sha, ...
    'expanded_content_sha256',content_sha, ...
    'roundtrip_verified',true, ...
    'independent_encodings_identical',true);
clear cleanup
local_cleanup({second,expanded});
end


function local_run_gzip(executable,source,destination)
command = sprintf('%s -n -9 -c -- %s > %s', ...
    local_shell_quote(executable),local_shell_quote(source), ...
    local_shell_quote(destination));
[status,output] = system(command);
if status ~= 0
    error('ver10_deterministic_gzip:CompressionFailed', ...
        'gzip failed with status %d: %s',status,strtrim(output));
end
end


function local_expand_gzip(executable,source,destination)
command = sprintf('%s -n -d -c -- %s > %s', ...
    local_shell_quote(executable),local_shell_quote(source), ...
    local_shell_quote(destination));
[status,output] = system(command);
if status ~= 0
    error('ver10_deterministic_gzip:ExpansionFailed', ...
        'gzip expansion failed with status %d: %s', ...
        status,strtrim(output));
end
end


function quoted = local_shell_quote(value)
value = char(value);
quoted = ['''',strrep(value,'''','''"''"'''),''''];
end


function local_cleanup(files)
for k = 1:numel(files)
    if exist(files{k},'file')
        delete(files{k});
    end
end
end
