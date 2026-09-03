function pool = ver10_ensure_local_pool(workers)
%VER10_ENSURE_LOCAL_POOL Start a local pool in an explicit job workspace.
%
% Set VER10_PARALLEL_JOB_STORAGE for production runs.  This prevents the
% local-cluster profile from silently writing worker metadata below a full
% home directory.  Interactive callers may omit it and retain MATLAB's
% normal local-profile behaviour.

if ~(isnumeric(workers) && isreal(workers) && isscalar(workers) ...
        && isfinite(workers) && workers >= 1 && workers == floor(workers))
    error('ver10_ensure_local_pool:BadWorkerCount', ...
        'workers must be a positive integer.');
end

storage = strtrim(getenv('VER10_PARALLEL_JOB_STORAGE'));
if ~isempty(storage)
    if ~local_is_absolute(storage)
        error('ver10_ensure_local_pool:RelativeStoragePath', ...
            'VER10_PARALLEL_JOB_STORAGE must be an absolute path.');
    end
    if ~isfolder(storage)
        [ok,message] = mkdir(storage);
        if ~ok
            error('ver10_ensure_local_pool:StorageCreateFailed', ...
                'Could not create parallel job storage: %s',message);
        end
    end
end

pool = gcp('nocreate');
if ~isempty(pool)
    reuse = pool.NumWorkers == workers;
    if reuse && ~isempty(storage)
        try
            reuse = strcmp(char(pool.Cluster.JobStorageLocation),storage);
        catch
            reuse = false;
        end
    end
    if reuse
        return;
    end
    delete(pool);
end

if isempty(storage)
    pool = parpool('local',workers);
else
    cluster = parcluster('local');
    cluster.JobStorageLocation = storage;
    pool = parpool(cluster,workers);
    if ~strcmp(char(pool.Cluster.JobStorageLocation),storage)
        delete(pool);
        error('ver10_ensure_local_pool:StorageMismatch', ...
            'The local pool did not use the requested job storage.');
    end
end
end


function tf = local_is_absolute(pathname)
tf = startsWith(pathname,filesep) ...
    || ~isempty(regexp(pathname,'^[A-Za-z]:[\\/]','once'));
end
