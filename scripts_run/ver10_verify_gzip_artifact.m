function certificate = ver10_verify_gzip_artifact( ...
    archive,expected_archive_sha256,expected_content_sha256)
%VER10_VERIFY_GZIP_ARTIFACT Verify a manifest-bound gzip and its contents.
%
% The archive is hashed before and after expansion.  Expansion goes to a
% fresh sibling temporary file, whose SHA-256 must equal the manifest's
% expanded-content digest.  This function never trusts gzip metadata or a
% previously computed hash.

archive = char(archive);
expected_archive_sha256 = char(expected_archive_sha256);
expected_content_sha256 = char(expected_content_sha256);
if exist(archive,'file') ~= 2 ...
        || isempty(regexp(expected_archive_sha256, ...
            '^[0-9a-f]{64}$','once')) ...
        || isempty(regexp(expected_content_sha256, ...
            '^[0-9a-f]{64}$','once'))
    error('ver10_verify_gzip_artifact:BadInput', ...
        'The gzip path or an expected SHA-256 is invalid.');
end
archive_sha_before = ver10_file_sha256(archive);
if ~strcmp(archive_sha_before,expected_archive_sha256)
    error('ver10_verify_gzip_artifact:ArchiveHashMismatch', ...
        'The gzip archive does not match its manifest SHA-256.');
end

directory = fileparts(archive);
if isempty(directory)
    directory = pwd;
end
expanded = [tempname(directory),'.expanded'];
cleanup = onCleanup(@() local_delete_if_present(expanded));
gzip_command = getenv('GZIP_COMMAND');
if isempty(gzip_command)
    gzip_command = 'gzip';
end
command = sprintf('%s -n -d -c -- %s > %s', ...
    local_shell_quote(gzip_command),local_shell_quote(archive), ...
    local_shell_quote(expanded));
[status,output] = system(command);
if status ~= 0
    error('ver10_verify_gzip_artifact:ExpansionFailed', ...
        'gzip expansion failed with status %d: %s', ...
        status,strtrim(output));
end
content_sha256 = ver10_file_sha256(expanded);
archive_sha_after = ver10_file_sha256(archive);
if ~strcmp(content_sha256,expected_content_sha256)
    error('ver10_verify_gzip_artifact:ContentHashMismatch', ...
        'The expanded gzip content has the wrong SHA-256.');
end
if ~strcmp(archive_sha_before,archive_sha_after) ...
        || ~strcmp(archive_sha_after,expected_archive_sha256)
    error('ver10_verify_gzip_artifact:ArchiveChangedDuringRead', ...
        'The gzip archive changed while it was being verified.');
end

certificate = struct( ...
    'schema','lowerboundsineq.gzip_artifact_verification.v1', ...
    'archive_sha256',archive_sha_after, ...
    'expanded_content_sha256',content_sha256, ...
    'archive_rehashed_after_expansion',true);
clear cleanup
local_delete_if_present(expanded);
end


function quoted = local_shell_quote(value)
value = char(value);
quoted = ['''',strrep(value,'''','''"''"'''),''''];
end


function local_delete_if_present(filename)
if exist(filename,'file')
    delete(filename);
end
end
