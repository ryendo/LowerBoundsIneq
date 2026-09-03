function test_intlab_init_lock_alias()
%TEST_INTLAB_INIT_LOCK_ALIAS Aliases of one tree must share one lock file.

if ~isunix
    fprintf('test_intlab_init_lock_alias: skipped outside POSIX.\n');
    return
end
root = tempname;
alias = [root,'-alias'];
mkdir(root);
cleanup = onCleanup(@() local_cleanup(root,alias));
[first_lock,first_path] = ver10_acquire_intlab_init_lock(root);
lock_cleanup = onCleanup(@() local_delete_file(first_path));
clear first_lock
[status,message] = system(sprintf('ln -s %s %s', ...
    local_shell_quote(root),local_shell_quote(alias)));
assert(status == 0,message);
[second_lock,second_path] = ver10_acquire_intlab_init_lock(alias);
clear second_lock
assert(strcmp(first_path,second_path), ...
    'A symlink alias selected a different initialization lock.');
clear lock_cleanup
local_delete_file(first_path);
clear cleanup
local_cleanup(root,alias);
fprintf('test_intlab_init_lock_alias: PASS (canonical shared lock path)\n');
end


function local_cleanup(root,alias)
system(sprintf('unlink %s >/dev/null 2>&1',local_shell_quote(alias)));
if isfolder(root)
    try
        rmdir(root,'s');
    catch
    end
end
end


function local_delete_file(filename)
if isfile(filename)
    try
        delete(filename);
    catch
    end
end
end


function quoted = local_shell_quote(value)
quoted = ['''',strrep(char(value),'''','''"''"'''),''''];
end
