function test_runtime_tree_hash_policy()
%TEST_RUNTIME_TREE_HASH_POLICY Exercise the exact MAT-v5 normalization gate.

root = tempname;
mkdir(root);
cleanup = onCleanup(@() local_remove_tree(root));
cache_name = 'Matlab_fixture_Intlab_Version_fixture.mat';
cache_path = fullfile(root,cache_name);
ordinary_path = fullfile(root,'ordinary.m');
local_write_binary(ordinary_path,uint8('proof-critical-source-v1'));
base_cache = local_mat_v5_fixture('Thu Sep  3 12:34:56 2026', ...
    'GLNXA64',uint8(0:63));
local_write_binary(cache_path,base_cache);

[tree1,count1,bytes1,records1] = ...
    ver10_directory_sha256(root,{cache_name});
assert(count1 == 2 && bytes1 == numel(base_cache)+24);
assert(numel(records1) == 1);
assert(strcmp(records1.relative_path,cache_name));
assert(records1.bytes == numel(base_cache));
assert(strcmp(records1.normalization, ...
    'mat-v5-created-on-24-character-timestamp-zeroed-v1'));
assert(~isempty(regexp(records1.normalized_sha256, ...
    '^[0-9a-f]{64}$','once')));

% Only the 24-character creation timestamp is volatile.
timestamp_only = local_mat_v5_fixture( ...
    'Fri Sep  4 23:45:01 2026','GLNXA64',uint8(0:63));
local_write_binary(cache_path,timestamp_only);
[tree2,count2,bytes2,records2] = ...
    ver10_directory_sha256(root,{cache_name});
assert(strcmp(tree1,tree2));
assert(count1 == count2 && bytes1 == bytes2);
assert(strcmp(records1.normalized_sha256, ...
    records2.normalized_sha256));

% Platform text, MAT control bytes, and serialized payload remain bound.
platform_changed = local_mat_v5_fixture( ...
    'Fri Sep  4 23:45:01 2026','GLNXB64',uint8(0:63));
local_write_binary(cache_path,platform_changed);
[platform_tree,~,~,platform_record] = ...
    ver10_directory_sha256(root,{cache_name});
assert(~strcmp(tree1,platform_tree));
assert(~strcmp(records1.normalized_sha256, ...
    platform_record.normalized_sha256));

payload_changed = timestamp_only;
payload_changed(end) = bitxor(payload_changed(end),uint8(1));
local_write_binary(cache_path,payload_changed);
[payload_tree,~,~,payload_record] = ...
    ver10_directory_sha256(root,{cache_name});
assert(~strcmp(tree1,payload_tree));
assert(~strcmp(records1.normalized_sha256, ...
    payload_record.normalized_sha256));

invalid_month = timestamp_only;
description = char(invalid_month(1:116));
month_start = strfind(description,'Sep');
assert(isscalar(month_start));
invalid_month(month_start:month_start+2) = uint8('Foo');
local_write_binary(cache_path,invalid_month);
local_assert_error(@() ver10_directory_sha256(root,{cache_name}), ...
    'ver10_directory_sha256:BadMatV5File');

invalid_endian = timestamp_only;
invalid_endian(127:128) = uint8('XX');
local_write_binary(cache_path,invalid_endian);
local_assert_error(@() ver10_directory_sha256(root,{cache_name}), ...
    'ver10_directory_sha256:BadMatV5File');

% Ordinary files are never normalized, even when the byte count is fixed.
local_write_binary(cache_path,timestamp_only);
local_write_binary(ordinary_path,uint8('proof-critical-source-v2'));
ordinary_tree = ver10_directory_sha256(root,{cache_name});
assert(~strcmp(tree1,ordinary_tree));

% Malformed MAT files and unsafe or non-text requested paths fail closed.
local_write_binary(cache_path,uint8(zeros(1,128)));
local_assert_error(@() ver10_directory_sha256(root,{cache_name}), ...
    'ver10_directory_sha256:BadMatV5File');
local_assert_error(@() ver10_directory_sha256(root,{'../escape.mat'}), ...
    'ver10_directory_sha256:BadNormalizedPath');
local_assert_error(@() ver10_directory_sha256(root,{17}), ...
    'ver10_directory_sha256:BadNormalizedFiles');

fprintf(['test_runtime_tree_hash_policy: PASS ', ...
    '(timestamp normalized; platform/payload/source tamper rejected)\n']);
clear cleanup
local_remove_tree(root);
end


function bytes = local_mat_v5_fixture(timestamp,platform,payload)
description = [ ...
    'MATLAB 5.0 MAT-file, Platform: ',char(platform), ...
    ', Created on: ',char(timestamp)];
assert(numel(timestamp) == 24 && numel(description) <= 116);
bytes = uint8(zeros(1,128+numel(payload)));
bytes(1:116) = uint8(' ');
bytes(1:numel(description)) = uint8(description);
bytes(125:126) = uint8([0,1]);
bytes(127:128) = uint8('IM');
bytes(129:end) = uint8(payload);
end


function local_write_binary(filename,bytes)
fid = fopen(filename,'wb');
assert(fid >= 0);
cleanup = onCleanup(@() fclose(fid));
written = fwrite(fid,uint8(bytes),'uint8');
assert(written == numel(bytes));
clear cleanup
end


function local_assert_error(callback,expected_identifier)
rejected = false;
try
    callback();
catch ME
    rejected = strcmp(ME.identifier,expected_identifier);
end
assert(rejected,['Expected rejection ',expected_identifier,'.']);
end


function local_remove_tree(directory)
if exist(directory,'dir')
    try
        rmdir(directory,'s');
    catch
    end
end
end
