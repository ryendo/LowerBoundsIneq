function my_intlab_config_worker()
%MY_INTLAB_CONFIG_WORKER  Per-parpool-worker INTLAB init for this project.
%   Called from parfevalOnAll to set up each worker's path and INTLAB state.
%   Unlike my_intlab_config, this does not attempt to regenerate the INTLAB
%   .mat; that is done once by the client. Also defensively scrubs other
%   Intlab installations that may appear on the shared system path.

    cur = strsplit(path, pathsep);
    bad = {};
    this_proj = fileparts(fileparts(mfilename('fullpath')));
    ok_intlab = fullfile(this_proj,'Intlab_V12');
    for k = 1:numel(cur)
        p = cur{k};
        if (contains(p, 'Intlab_V12') || contains(p, 'INTLAB')) && ...
           ~startsWith(p, ok_intlab)
            bad{end+1} = p; %#ok<AGROW>
        end
    end
    for k = 1:numel(bad), try, rmpath(bad{k}); catch, end, end

    addpath(this_proj);
    my_intlab_config();
end
