function clusters = auto_cluster_eigenvalues(lam, relative_threshold)
% AUTO_CLUSTER_BY_WIDTH: Automatically groups eigenvalues into clusters
% based on the total width of the cluster.
%
% A cluster grows as long as its width (λ_max - λ_min) does not exceed a
% threshold defined relative to its starting eigenvalue.
%
% INPUTS:
%   lam:                (Interval Vector) Rigorous bounds for eigenvalues.
%   relative_threshold: (Scalar) The maximum allowed relative width of a cluster.
%                       Example: 0.02 means a cluster's width cannot exceed 2%
%                       of its starting eigenvalue's value.
%
% OUTPUT:
%   clusters:           (Cell Array) The determined eigenvalue clusters.

num_eigs = length(lam);
clusters = {};
start_idx = 1;

while start_idx <= num_eigs
    % Define the absolute width threshold for this new cluster based on its
    % starting eigenvalue. We use the upper bound for a robust criterion.
    threshold = I_sup(lam(start_idx)) * relative_threshold;
    
    end_idx = start_idx;
    for j = (start_idx + 1):num_eigs
        % Calculate the width of the potential cluster from start_idx to j
        cluster_width = I_sup(lam(j)) - I_inf(lam(start_idx));
        
        if cluster_width > threshold
            % The width exceeds the threshold, so the cluster ends at j-1.
            break;
        end
        % If width is within the threshold, expand the cluster
        end_idx = j;
    end
    
    % Finalize and store the current cluster
    clusters{end+1} = start_idx:end_idx;
    
    % Move to the start of the next cluster
    start_idx = end_idx + 1;
end

% fprintf('Automatic clustering by width determined %d clusters.\n', length(clusters));
end