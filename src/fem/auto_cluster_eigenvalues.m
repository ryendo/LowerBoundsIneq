function clusters = auto_cluster_eigenvalues(lam, relative_gap_threshold)
%AUTO_CLUSTER_EIGENVALUES  Group eigenvalue enclosures into eigenspaces.
%
% A cluster grows while consecutive eigenvalue enclosures are not provably
% separated.  This is the safe choice for eigenspace estimates: if two
% enclosures overlap, a one-dimensional eigenfunction error is not a
% well-posed certified quantity, and the whole overlapping block must be
% treated as one eigenspace.
%
% relative_gap_threshold optionally merges very small positive gaps as a
% conservative near-degeneracy guard.  Use 0 to merge only overlapping
% intervals.  Existing callers use 0.01.
%
% INPUTS:
%   lam:                (Interval Vector) Rigorous bounds for eigenvalues.
%   relative_gap_threshold:
%                       (Scalar) Merge if the positive gap before the next
%                       enclosure is no larger than this fraction of the
%                       cluster scale.  Example: 0.01 means 1%.
%
% OUTPUT:
%   clusters:           (Cell Array) The determined eigenvalue clusters.

num_eigs = length(lam);
clusters = {};
start_idx = 1;

while start_idx <= num_eigs
    end_idx = start_idx;
    cluster_sup = I_sup(lam(start_idx));
    cluster_scale = max(1, abs(I_mid(lam(start_idx))));

    for j = (start_idx + 1):num_eigs
        gap = I_inf(lam(j)) - cluster_sup;
        gap_tolerance = relative_gap_threshold * cluster_scale;

        if gap > gap_tolerance
            break;
        end

        end_idx = j;
        cluster_sup = max(cluster_sup, I_sup(lam(j)));
        cluster_scale = max(cluster_scale, abs(I_mid(lam(j))));
    end

    clusters{end+1} = start_idx:end_idx;
    start_idx = end_idx + 1;
end

end
