function [Est_grad, delta_a_sq, delta_b_sq, info] = calc_grad_error_bounds( ...
    lam, lam_h, uh_list, A, B, clusters)

    % normalize clusters to numeric row vectors
    num_clusters = numel(clusters);
    cl = cell(num_clusters,1);
    for k = 1:num_clusters
        idx = clusters{k};
        if iscell(idx), idx = idx{1}; end
        cl{k} = idx(:).';
    end
    clusters = cl;

    num_eigs = length(lam);

    delta_a_sq = I_intval(I_ones(num_clusters,1));
    delta_b_sq = I_intval(I_ones(num_clusters,1));

    info = struct();
    info.ok     = false(num_clusters,1);
    info.reason = repmat({''}, num_clusters,1);
    info.rho    = I_intval(NaN(num_clusters,1));

    for k = 1:num_clusters
        idx_k = clusters{k};
        n_k   = idx_k(1);
        N_k   = idx_k(end);

        next_index = N_k + 1;
        if next_index > num_eigs
            info.reason{k} = 'No λ_{N_k+1} available (last cluster). Lemma 1 not applicable -> Inf.';
            continue;
        end

        rho = lam(next_index);
        info.rho(k) = rho;

        % Theorem requires λ_{n_k} < rho <= λ_{N_k+1}.
        if ~(lam(n_k) < rho)
            info.reason{k} = 'Cannot certify λ_{n_k} < rho via interval comparison -> Inf.';
            continue;
        end

        % Need denominators positive: rho - λ_{n_k} > 0 and λ_{N_k,h} > 0
        den_gap = rho - lam(n_k);
        if ~(I_intval(0) < den_gap)
            info.reason{k} = 'Cannot certify (rho - λ_{n_k}) > 0 via interval comparison -> Inf.';
            continue;
        end

        % λ_{N_k,h} := max over E_k^h Rayleigh quotient.
        lamNkh = max(I_sup(lam_h(idx_k)));
        if ~(I_intval(0) < lamNkh)
            info.reason{k} = 'Cannot certify λ_{N_k,h} > 0 -> Inf.';
            continue;
        end

        % --- theta term
        theta_a = I_intval(0);
        theta_b = I_intval(0);

        for l = 1:k-1
            idx_l = clusters{l};
            n_l   = idx_l(1);

            % If previous cluster is not certified, theta becomes Inf
            if isinf(delta_a_sq(l)) || isinf(delta_b_sq(l))
                theta_a = I_intval(Inf);
                theta_b = I_intval(Inf);
                break;
            end

            eps_a = epsilon_subspaces_bound(uh_list(:,idx_l), uh_list(:,idx_k), A);
            eps_b = epsilon_subspaces_bound(uh_list(:,idx_l), uh_list(:,idx_k), B);

            % Need λ_{n_l} > 0 to define (rho-λ)/λ
            if ~(I_intval(0) < lam(n_l))
                theta_a = I_intval(Inf);
                theta_b = I_intval(Inf);
                break;
            end

            factor_a = (rho - lam(n_l)) ./ lam(n_l);
            factor_b = (rho - lam(n_l));

            da_l = sqrt(delta_a_sq(l));
            db_l = sqrt(delta_b_sq(l));

            theta_a = theta_a + factor_a .* (eps_a + da_l).^2;
            theta_b = theta_b + factor_b .* (eps_b + db_l).^2;
        end

        if isinf(theta_a) || isinf(theta_b)
            info.reason{k} = 'theta became Inf (previous uncertified cluster or λ_{n_l} not positive).';
            continue;
        end

        % --- Liu–Vejchodský Lemma 1 (all intval, no clipping)
        lam_nk_inf = I_inf(lam(n_k));
        lam_nk_sup = I_sup(lam(n_k));

        den_a = lamNkh .* (rho - lam_nk_sup);
        den_b = (rho - lam_nk_sup);

        if ~(I_intval(0) < den_a) || ~(I_intval(0) < den_b)
            info.reason{k} = 'Cannot certify denominators > 0 -> Inf.';
            continue;
        end

        num_a = rho.*(lamNkh - lam_nk_inf) + lam_nk_sup.*lamNkh.*theta_a;
        num_b = (lamNkh - lam_nk_inf) + theta_b;

        delta_a_sq(k) = num_a ./ den_a;
        delta_b_sq(k) = num_b ./ den_b;

        info.ok(k) = true;
        info.reason{k} = '';
    end

    % --- Assign per-eigenfunction cluster-wise (subspace bound)
    Est_grad = I_intval(Inf(num_eigs,1));
    for k = 1:num_clusters
        idx_k = clusters{k};
        if ~isinf(delta_a_sq(k))
            Est_grad(idx_k) = sqrt(delta_a_sq(k));
        else
            Est_grad(idx_k) = I_intval(Inf);
        end
    end
end


% =========================================================================
% epsilon_subspaces_bound:
% Lemma-2 style bound without any orthonormalization:
%
% Let G = U^T M U, H = V^T M V, F = U^T M V.
% If ||I-G||_2 <= eta_G < 1 and ||I-H||_2 <= eta_H < 1 and ||F^T F||_2 <= eta_F,
% then eps^2 <= eta_F / ((1-eta_G)(1-eta_H)).
%
% If conditions are not certifiable, return [0,1].
% =========================================================================
function eps = epsilon_subspaces_bound(U, V, M)
    if isempty(U) || isempty(V)
        eps = I_infsup(0,0);
        return;
    end

    if ~isa(U,'intval'), U = I_intval(U); end
    if ~isa(V,'intval'), V = I_intval(V); end
    if ~isa(M,'intval'), M = I_intval(M); end

    % --- Orthonormalize bases w.r.t. M (THIS is the key fix)
    [Uo, okU] = orthonormalize_wrt_M(U, M);
    [Vo, okV] = orthonormalize_wrt_M(V, M);
    if ~(okU && okV)
        eps = I_infsup(0,1);
        return;
    end

    % --- Rebuild Gram blocks after orthonormalization
    G = Uo' * M * Uo;
    H = Vo' * M * Vo;
    F = Uo' * M * Vo;

    m = size(G,1);
    p = size(H,1);

    eta_G = norm(I_intval(eye(m)) - G, 2);
    eta_H = norm(I_intval(eye(p)) - H, 2);

    % Use ||F^T F||_2 <= ||F||_2^2 (less blow-up than forming F'*F)
    eta_F = norm(F, 2)^2;

    if ~(eta_G < 1) || ~(eta_H < 1)
        eps = I_infsup(0,1);
        return;
    end

    eps2 = eta_F ./ ((1 - eta_G).*(1 - eta_H));
    eps  = sqrt(eps2);
end



% =========================================================================
% Orthonormalize columns of U w.r.t. M
% =========================================================================
function [Uo, ok] = orthonormalize_wrt_M(U, M)
    ok = true;

    % Interval Gram matrix (rigorous)
    Gint = U' * M * U;                 % intval (m x m)

    % Rigorous symmetrization in interval arithmetic
    Gsym = (Gint + Gint')/2;

    % Use midpoint only to build a preconditioner (R)
    Gmid = mid(Gsym);                  % double symmetric

    % Cholesky on the midpoint (preconditioner)
    [R, pflag] = chol(Gmid);
    if pflag ~= 0
        ok = false;
        Uo = U;
        return;
    end

    % Transform basis
    Uo = U / R;
end
