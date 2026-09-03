function [cleanup,lock_path] = ver10_acquire_intlab_init_lock(intlab_root)
%VER10_ACQUIRE_INTLAB_INIT_LOCK Serialize writes to INTLAB's data MAT.
%
% STARTINTLAB rewrites a shared, proof-critical platform-data MAT file on
% every invocation.  Parallel MATLAB workers must therefore not initialize
% INTLAB concurrently.  A Java NIO advisory lock is used because it is
% released by the operating system if a MATLAB process is killed.

intlab_root = char(intlab_root);
if isempty(intlab_root) || ~isfolder(intlab_root)
    error('ver10_acquire_intlab_init_lock:BadIntlabRoot', ...
        'INTLAB_ROOT must name an existing directory.');
end
try
    root_file = javaObject('java.io.File',intlab_root);
    intlab_root = char(root_file.getCanonicalPath());
catch ME
    error('ver10_acquire_intlab_init_lock:CannotCanonicalizeRoot', ...
        'Cannot canonicalize INTLAB_ROOT: %s',ME.message);
end
[parent,name] = fileparts(regexprep(intlab_root,[regexptranslate( ...
    'escape',filesep),'+$'],''));
if isempty(parent) || isempty(name)
    error('ver10_acquire_intlab_init_lock:BadIntlabRoot', ...
        'Cannot derive a safe INTLAB initialization lock path.');
end
lock_path = fullfile(parent,['.',name,'.ver10-startintlab.lock']);

try
    random_access = javaObject('java.io.RandomAccessFile',lock_path,'rw');
    channel = random_access.getChannel();
catch ME
    error('ver10_acquire_intlab_init_lock:CannotOpenLock', ...
        'Cannot open INTLAB initialization lock: %s',ME.message);
end

lock = [];
started = tic;
while isempty(lock)
    try
        lock = channel.tryLock();
    catch ME
        if ~contains(char(ME.message),'OverlappingFileLockException')
            local_close_unlocked(channel,random_access);
            error('ver10_acquire_intlab_init_lock:CannotAcquireLock', ...
                'Cannot acquire INTLAB initialization lock: %s', ...
                ME.message);
        end
    end
    if ~isempty(lock)
        break
    end
    if toc(started) > 300
        local_close_unlocked(channel,random_access);
        error('ver10_acquire_intlab_init_lock:Timeout', ...
            'Timed out waiting for the INTLAB initialization lock.');
    end
    pause(0.05);
end
cleanup = onCleanup(@() local_release(lock,channel,random_access));
end


function local_release(lock,channel,random_access)
try
    if ~isempty(lock) && lock.isValid()
        lock.release();
    end
catch
end
local_close_unlocked(channel,random_access);
end


function local_close_unlocked(channel,random_access)
try
    channel.close();
catch
end
try
    random_access.close();
catch
end
end
