function source = ver10_assert_source_unchanged(project_root,source)
%VER10_ASSERT_SOURCE_UNCHANGED Fail if a run's Git source changed in flight.
%
% The returned structure is the start-of-run provenance enriched with the
% independent end-of-run check.  Production artifacts can therefore state
% both the commit used to form checkpoint fingerprints and that the same
% worktree state was present when the certificate was assembled.

finish = ver10_source_provenance(project_root);
if source.git_status_available ~= finish.git_status_available
    error('ver10_assert_source_unchanged:SourceChangedDuringRun', ...
        'Git provenance availability changed during the run.');
end
if source.git_status_available
    same_state = strcmp(source.git_commit,finish.git_commit) ...
        && source.git_dirty == finish.git_dirty ...
        && strcmp(source.git_status_porcelain,finish.git_status_porcelain);
    if ~same_state
        error('ver10_assert_source_unchanged:SourceChangedDuringRun', ...
            ['The Git commit or worktree state changed during the run.  ', ...
             'The resulting checkpoints cannot be published as one-source ', ...
             'certificate.']);
    end
end

source.provenance_rechecked_after_run = true;
source.git_commit_after_run = finish.git_commit;
source.git_dirty_after_run = finish.git_dirty;
source.source_unchanged_after_run = source.git_status_available;
end
