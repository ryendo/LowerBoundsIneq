function [digest,file_count,total_bytes,normalized_records] = ...
    ver10_directory_sha256(directory,normalized_mat_relative_names, ...
        normalization_lock_factory)
%VER10_DIRECTORY_SHA256 Hash a tree with an explicit MAT-v5 header policy.
%
% INTLAB rewrites its platform data MAT file at every initialization even
% when the proof-critical INTLAB_CONST payload is unchanged.  MATLAB's
% level-5 MAT writer embeds a creation timestamp in the 116-byte
% descriptive header.  For the explicitly named files only, this routine
% replaces exactly the 24-character value following "Created on: " by
% zeroes before hashing.  The rest of the description (including platform),
% path, byte count, remaining MAT header, and complete serialized payload
% remain part of both the tree digest and the per-file normalized digest.

% No filename pattern is interpreted here: callers must supply exact
% relative paths.  This keeps the normalization boundary auditable and
% prevents an unrelated MAT file from being silently normalized.

if nargin < 2
    normalized_mat_relative_names = {};
end
if nargin < 3
    normalization_lock_factory = [];
end
if ~isempty(normalization_lock_factory) ...
        && ~isa(normalization_lock_factory,'function_handle')
    error('ver10_directory_sha256:BadLockFactory', ...
        'The normalization lock factory must be a function handle.');
end
directory = char(directory);
if isempty(directory) || ~isfolder(directory)
    digest = '';
    file_count = 0;
    total_bytes = 0;
    normalized_records = local_empty_records();
    return
end

requested = local_validate_relative_names( ...
    normalized_mat_relative_names);
entries = dir(fullfile(directory,'**','*'));
entries = entries(~[entries.isdir]);
if isempty(entries)
    if ~isempty(requested)
        error('ver10_directory_sha256:MissingNormalizedFile', ...
            'A requested normalized MAT file is absent from the tree.');
    end
    digest = local_empty_tree_digest();
    file_count = 0;
    total_bytes = 0;
    normalized_records = local_empty_records();
    return
end

full_names = arrayfun(@(entry) fullfile(entry.folder,entry.name), ...
    entries,'UniformOutput',false);
root_prefix = [char(directory),filesep];
relative_names = full_names;
for k = 1:numel(relative_names)
    local_name = full_names{k};
    if ~startsWith(local_name,root_prefix)
        error('ver10_directory_sha256:BadTreeEntry', ...
            'Tree entry is outside the recorded root: %s',local_name);
    end
    relative_names{k} = strrep( ...
        local_name(numel(root_prefix)+1:end),filesep,'/');
end
[relative_names,order] = sort(relative_names);
full_names = full_names(order);
entries = entries(order);

missing = setdiff(requested,relative_names);
if ~isempty(missing)
    error('ver10_directory_sha256:MissingNormalizedFile', ...
        'Requested normalized MAT file is absent: %s',missing{1});
end

tree_engine = local_new_engine();
local_digest_update(tree_engine,uint8( ...
    'ver10-directory-sha256-v2-mat5-created-on-normalization'));
local_digest_update(tree_engine,uint8(0));
for k = 1:numel(requested)
    local_digest_update(tree_engine,uint8('normalized-mat-v5'));
    local_digest_update(tree_engine,uint8(0));
    local_digest_update(tree_engine,unicode2native(requested{k},'UTF-8'));
    local_digest_update(tree_engine,uint8(0));
end

normalized_records = local_empty_records();
for k = 1:numel(full_names)
    relative_name = relative_names{k};
    size_value = double(entries(k).bytes);
    normalize_header = ismember(relative_name,requested);

    local_digest_update(tree_engine, ...
        unicode2native(relative_name,'UTF-8'));
    local_digest_update(tree_engine,uint8(0));
    local_digest_update(tree_engine,uint8(sprintf('%d',size_value)));
    local_digest_update(tree_engine,uint8(0));

    file_engine = [];
    if normalize_header
        if size_value < 128
            error('ver10_directory_sha256:BadMatV5File', ...
                'Normalized MAT-v5 file is shorter than 128 bytes: %s', ...
                relative_name);
        end
        file_engine = local_new_engine();
    end

    normalization_lock = [];
    if normalize_header && ~isempty(normalization_lock_factory)
        normalization_lock = normalization_lock_factory();
    end
    fid = fopen(full_names{k},'rb');
    if fid < 0
        clear normalization_lock
        error('ver10_directory_sha256:CannotOpenFile', ...
            'Cannot open tree file %s.',full_names{k});
    end
    cleanup = onCleanup(@() fclose(fid));
    offset = 0;
    first_chunk = true;
    while true
        bytes = fread(fid,1024*1024,'*uint8');
        if isempty(bytes)
            break
        end
        if normalize_header
            if first_chunk
                bytes = local_normalize_mat_v5_header( ...
                    bytes,relative_name);
                first_chunk = false;
            end
            local_digest_update(file_engine,bytes);
        end
        local_digest_update(tree_engine,bytes);
        offset = offset+numel(bytes);
    end
    clear cleanup
    clear normalization_lock
    local_digest_update(tree_engine,uint8(0));

    if offset ~= size_value
        error('ver10_directory_sha256:FileChangedWhileHashing', ...
            'Tree file changed size while it was hashed: %s',relative_name);
    end
    if normalize_header
        record = struct( ...
            'relative_path',relative_name, ...
            'bytes',size_value, ...
            'normalization', ...
                'mat-v5-created-on-24-character-timestamp-zeroed-v1', ...
            'normalized_sha256',local_hex_digest(file_engine));
        normalized_records(end+1) = record; %#ok<AGROW>
    end
end

digest = local_hex_digest(tree_engine);
file_count = numel(entries);
total_bytes = sum(double([entries.bytes]));
end


function requested = local_validate_relative_names(value)
if isempty(value)
    requested = {};
elseif ischar(value) || (isstring(value) && isscalar(value))
    requested = {char(value)};
elseif iscell(value)
    requested = cell(size(value));
    for k = 1:numel(value)
        item = value{k};
        if ~(ischar(item) || (isstring(item) && isscalar(item)))
            error('ver10_directory_sha256:BadNormalizedFiles', ...
                'Every normalized MAT path must be a text scalar.');
        end
        requested{k} = char(item);
    end
elseif isstring(value)
    requested = cellstr(value(:));
else
    error('ver10_directory_sha256:BadNormalizedFiles', ...
        'Normalized MAT paths must be text values.');
end
for k = 1:numel(requested)
    requested{k} = strrep(requested{k},'\','/');
    parts = strsplit(requested{k},'/');
    if isempty(requested{k}) || startsWith(requested{k},'/') ...
            || ~isempty(regexp(requested{k},'^[A-Za-z]:','once')) ...
            || any(strcmp(parts,'.')) || any(strcmp(parts,'..'))
        error('ver10_directory_sha256:BadNormalizedPath', ...
            'A normalized MAT path is not a safe relative path.');
    end
end
requested = sort(unique(requested));
end


function bytes = local_normalize_mat_v5_header(bytes,relative_name)
if numel(bytes) < 128 ...
        || ~strcmp(char(reshape(bytes(1:19),1,[])), ...
            'MATLAB 5.0 MAT-file')
    error('ver10_directory_sha256:BadMatV5File', ...
        ['Only the creation timestamp of a MATLAB level-5 MAT ', ...
         'file may be normalized: %s'],relative_name);
end
endian = char(reshape(bytes(127:128),1,[]));
version_bytes = reshape(bytes(125:126),1,[]);
valid_version = (strcmp(endian,'IM') ...
        && isequal(version_bytes,uint8([0,1]))) ...
    || (strcmp(endian,'MI') ...
        && isequal(version_bytes,uint8([1,0])));
if ~valid_version
    error('ver10_directory_sha256:BadMatV5File', ...
        'MAT-v5 version or endian bytes are invalid: %s',relative_name);
end
description = char(reshape(bytes(1:116),1,[]));
pattern = [ ...
    '^MATLAB 5\.0 MAT-file, Platform: [A-Za-z0-9_.-]+, Created on: ', ...
    '(Sun|Mon|Tue|Wed|Thu|Fri|Sat) ', ...
    '(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec) ', ...
    '( [1-9]|[12][0-9]|3[01]) ', ...
    '([01][0-9]|2[0-3]):[0-5][0-9]:[0-5][0-9] [0-9]{4} *$'];
if isempty(regexp(description,pattern,'once'))
    error('ver10_directory_sha256:BadMatV5File', ...
        'MAT-v5 description is not in the canonical MATLAB format: %s', ...
        relative_name);
end
marker = strfind(description,'Created on: ');
if numel(marker) ~= 1
    error('ver10_directory_sha256:BadMatV5File', ...
        'MAT-v5 creation marker is absent or ambiguous: %s', ...
        relative_name);
end
timestamp_start = marker+numel('Created on: ');
timestamp_end = timestamp_start+23;
bytes(timestamp_start:timestamp_end) = uint8(0);
end


function records = local_empty_records()
records = struct('relative_path',{},'bytes',{}, ...
    'normalization',{},'normalized_sha256',{});
end


function digest = local_empty_tree_digest()
engine = local_new_engine();
local_digest_update(engine,uint8( ...
    'ver10-directory-sha256-v2-mat5-created-on-normalization'));
local_digest_update(engine,uint8(0));
digest = local_hex_digest(engine);
end


function engine = local_new_engine()
try
    engine = javaMethod( ...
        'getInstance','java.security.MessageDigest','SHA-256');
catch ME
    error('ver10_directory_sha256:NoSHA256', ...
        'Cannot initialize the SHA-256 engine: %s',ME.message);
end
end


function local_digest_update(engine,bytes)
engine.update(typecast(uint8(bytes(:)),'int8'));
end


function digest = local_hex_digest(engine)
raw = typecast(engine.digest(),'uint8');
digest = lower(reshape(dec2hex(raw,2).',1,[]));
end
