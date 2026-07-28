function source = ver10_source_provenance(project_root)
%VER10_SOURCE_PROVENANCE Record the exact source tree used by a run.

old_directory = pwd;
cleanup = onCleanup(@() cd(old_directory));
cd(project_root);

[commit_status,commit_text] = system('git rev-parse HEAD');
[state_status,state_text] = system( ...
    'git status --porcelain --untracked-files=normal');

source = struct();
if commit_status == 0
    source.git_commit = strtrim(commit_text);
else
    source.git_commit = '';
end
source.git_status_available = commit_status == 0 && state_status == 0;
source.git_dirty = state_status ~= 0 || ~isempty(strtrim(state_text));
source.git_status_porcelain = strtrim(state_text);
source.provenance_captured_before_run = true;

clear cleanup
end
