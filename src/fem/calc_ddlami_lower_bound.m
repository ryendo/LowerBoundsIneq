function [lami, dlami, ddlam_i_lower_bound] = calc_ddlami_lower_bound(i, base_triangle, triangle, e_direction, N_spectral, N_LG, N_rho, fem_ord_LG)
% CALC_DDLAMI_LOWER_BOUND
% Computes a rigorous lower bound for the second-order directional shape
% derivative of the i-th Dirichlet eigenvalue on a triangle, using validated
% eigenvalue bounds and a FEM-based truncation of the spectral expansion.
%
% This implementation is aligned with the paper's spectral representation
% (Lemma "Spectral representation of the second-order shape derivative") and
% the computable certified lower bound (Theorem \ref{thm:main-theorem-est}):
%
%   ddot(lambda_i^p)
%   >=  widehat(ddot(lambda_{i,N}^p))
%       - ||ddot P||_2 ( sqrt(lambda_i^p) + sqrt(hat(lambda)_i^p) ) Est_a(u_i^p, hat(u)_i^p)
%       - 2 * sum_{k=1..N, k!=i}  A_{ik} / |lambda_k^p - lambda_i^p|
%
% where
%   A_{ik}
%   := ||dot P||_2^2
%      * ( Est_a(u_i^p,hat(u)_i^p) * sqrt(hat(lambda)_k^p) + sqrt(lambda_i^p) * Est_a(u_k^p,hat(u)_k^p) )
%      * ( sqrt(lambda_i^p*lambda_k^p) + sqrt(hat(lambda)_i^p*hat(lambda)_k^p) ).
%
% IMPORTANT:
% - The theorem assumes i < N (here N = N_spectral).
% - The bounds require that the eigenvalue enclosures for lambda_i and lambda_k
%   are separated so that (lambda_k - lambda_i) does not contain 0.

% fprintf('--- Start calc_ddlami_lower_bound for eigenvalue index i = %d ---\n', i);

% Theorem \ref{thm:main-theorem-est} requires i < N.
if i >= N_spectral
    error('The paper''s bound requires i < N_spectral. Please use N_spectral >= i+1.');
end

% We typically need at least N_spectral+1 eigenpairs; the extra "+2" is kept
% to match the original workflow (e.g., clustering / separation checks).
num_eigs_needed = N_spectral + 1;

% fprintf('  > Computing eigenvalue bounds and eigenfunctions (num_eigs_needed = %d)...\n', num_eigs_needed);

% lams_raw: rigorous interval enclosures of true eigenvalues on base_triangle.
% uh_list : corresponding FEM eigenvectors (used as discrete eigenfunctions).
% A,B     : stiffness/mass matrices, and A_xx,A_xy,A_yy for derivative bilinear forms.

format long infsup

[lams_raw, uh_list, A, B, A_xx, A_xy, A_yy] = ...
    calc_eigen_bounds_any_order_1k_wh(num_eigs_needed, base_triangle, N_LG, N_rho, fem_ord_LG);

% fprintf('    done.\n');

% Compute discrete eigenvalues (Rayleigh quotients) on base_triangle.
lams_h_raw = I_intval(zeros(num_eigs_needed, 1));
for k = 1:num_eigs_needed
    lams_h_raw(k) = (uh_list(:,k)' * A * uh_list(:,k)) / (uh_list(:,k)' * B * uh_list(:,k));
end

% --- Coordinate transformation and eigenvalue perturbation bounds (paper's Lemma on eigenvalue perturbation) ---
% fprintf('  > Applying coordinate transformation and perturbation factors...\n');

x_base   = base_triangle(5);
y_base   = base_triangle(6);
x_target = triangle(5);
y_target = triangle(6);

% Linear map S_{base,target} sends triangle^{base} to triangle^{target}.
% Here S_inv = S^{-1}. For the eigenvalue perturbation lemma we need S^{-1}S^{-T}.
S_inv     = [1, (x_base - x_target) / y_target; 0, y_base / y_target];
SinvSinvt = S_inv * S_inv';

% Bounds for eigenvalues of SinvSinvt (min/max), used as a common scaling interval.
scale_interval = norm(SinvSinvt);

% Rigorous bounds for true eigenvalues and discrete eigenvalues on the target triangle.
lams_interval  = lams_raw   * scale_interval;   % enclosures for lambda_k^p
lams_h         = lams_h_raw * scale_interval;   % proxy for hat(lambda)_k^p

% fprintf('    done.\n');

% --- Error estimators for eigenfunctions (Est_a) ---
% fprintf('  > Calculating error bounds (Est_grad) and clustering eigenvalues...\n');

clusters = auto_cluster_eigenvalues(lams_interval, 0.01);

% Est_grad(j) is used as Est_a(u_j^p, hat(u)_j^p) in the theorem.
Est_grad = calc_grad_error_bounds(lams_interval, lams_h, uh_list, A, B, clusters);

% fprintf('    done.\n');

% --- Shape derivative matrices dot P and ddot P (paper's Lemma 1st/2nd derivative formulas) ---
a = e_direction(1);
b = e_direction(2);

mat_dP   = [[0,               -a/y_target];         [-a/y_target,       -2*b/y_target]];
mat_ddP  = [[2*a^2/y_target^2, 4*a*b/y_target^2]; [4*a*b/y_target^2, 6*b^2/y_target^2]];

% Spectral (operator) norms as in the theorem.
norm_dP  = norm(mat_dP, 2);
norm_ddP = norm(mat_ddP, 2);

% --- Eigenvalue enclosure to be returned (use the transformed/target enclosure) ---
lami   = lams_interval(i);   % rigorous enclosure for lambda_i^p (target triangle)
lamih  = lams_h(i);          % proxy enclosure for hat(lambda)_i^p

%======================================================================
% STEP 1: FIRST-ORDER DIRECTIONAL DERIVATIVE (Hadamard formula + error bound)
%======================================================================
% fprintf('  > Computing first-order directional derivative...\n');

% Discrete approximation: (dot P ∇hat(u_i), ∇hat(u_i))_{triangle^p}.
dlamih = uh_list(:,i)' * ( ...
    mat_dP(1,1)*A_xx + (mat_dP(1,2)+mat_dP(2,1))*A_xy + mat_dP(2,2)*A_yy ) * uh_list(:,i);

% Certified a priori bound (same structure as the proof for the 2nd derivative term):
% |dot(lambda_i^p) - dot(lambda_{h,i}^p)|
% <= ||dot P||_2 * (sqrt(lambda_i^p) + sqrt(hat(lambda)_i^p)) * Est_a(u_i^p, hat(u)_i^p).
e_dlami = norm_dP * (sqrt(lami) + sqrt(lamih)) * Est_grad(i);

dlami = I_hull(dlamih - e_dlami, dlamih + e_dlami);

% fprintf('    done.\n');

%======================================================================
% STEP 2: COMPUTE widehat(ddot(lambda_{i,N}^p))
%======================================================================
% fprintf('  > Step 2: Computing widehat(ddot(lambda_{%d,%d}^p))...\n', i, N_spectral);

ddlam_i_N_h = I_intval('0');

% T1_h := (ddot P ∇hat(u_i), ∇hat(u_i))_{triangle^p}.
T1_h = uh_list(:,i)' * ( ...
    mat_ddP(1,1)*A_xx + (mat_ddP(1,2)+mat_ddP(2,1))*A_xy + mat_ddP(2,2)*A_yy ) * uh_list(:,i);

ddlam_i_N_h = ddlam_i_N_h + T1_h;

% fprintf('  > Adding truncated spectral sum (k=1..%d, k!=i)...\n', N_spectral);

for k = 1:N_spectral
    if k == i
        continue;
    end

    % num := (dot P ∇hat(u_i), ∇hat(u_k))_{triangle^p}.
    num = uh_list(:,i)' * ( ...
        mat_dP(1,1)*A_xx + (mat_dP(1,2)+mat_dP(2,1))*A_xy + mat_dP(2,2)*A_yy ) * uh_list(:,k);

    % den := (lambda_k^p - lambda_i^p), rigorously enclosed.
    den = lams_interval(k) - lams_interval(i);

    % For certification we must avoid division by an interval containing 0.
    if (I_inf(den) <= 0) && (I_sup(den) >= 0)
        error('Eigenvalue enclosures overlap: (lambda_%d - lambda_%d) contains 0. Cannot certify the spectral term.', k, i);
    end

    ddlam_i_N_h = ddlam_i_N_h + 2 * (abs(num)^2) / den;
end

% fprintf('    done.\n');

%======================================================================
% STEP 3: COMPUTE the theorem's error term and assemble the rigorous lower bound
%======================================================================
% fprintf('  > Step 3: Assembling the theorem-based rigorous lower bound...\n');

% E_approx is the full subtraction term in Theorem \ref{thm:main-theorem-est}:
%   ||ddot P||_2 (sqrt(lambda_i)+sqrt(hat lambda_i)) Est_a(u_i,hat u_i)
% + 2 sum_{k!=i} A_ik / |lambda_k-lambda_i|.
E_approx = I_intval('0');

% First (ddot P)-term in the theorem.
E_approx = E_approx + norm_ddP * (sqrt(lami) + sqrt(lamih)) * Est_grad(i);

% Spectral-sum error terms in the theorem.
for k = 1:N_spectral
    if k == i
        continue;
    end

    den = lams_interval(k) - lams_interval(i);

    % Same separation check as above (also needed for |den| in the theorem).
    if (I_inf(den) <= 0) && (I_sup(den) >= 0)
        error('Eigenvalue enclosures overlap: (lambda_%d - lambda_%d) contains 0. Cannot certify the error term.', k, i);
    end

    % Theorem's A_ik:
    A_ik = norm_dP^2 ...
        * ( Est_grad(i)*sqrt(lams_h(k)) + sqrt(lami)*Est_grad(k) ) ...
        * ( sqrt(lami*lams_interval(k)) + sqrt(lamih*lams_h(k)) );

    % Theorem uses 2 * A_ik / |lambda_k - lambda_i|.
    error_k = 2 * (A_ik / abs(den));

    E_approx = E_approx + error_k;
end

% Final rigorous lower bound
ddlam_i_lower_bound = ddlam_i_N_h - E_approx;


% fprintf('    done.\n');

try
    lamNp1   = lams_interval(N_spectral + 1);
    den_tail = lamNp1 - lams_interval(i);   % (lambda_{N+1}^p - lambda_i^p)

    if (I_inf(den_tail) <= 0)
        warning('calc_ddlami_lower_bound:TailBoundSkipped', ...
            '%s', 'Tail bound skipped: lambda_{N+1}-lambda_i is not provably positive.');
    else
        norm_dP_F = norm(mat_dP, 'fro');    % ||dot P||_F
        R_tail_ub = 2 * (norm_dP_F^2) * (lami^2) / den_tail;        

        if I_inf(R_tail_ub) < 0
            warning('calc_ddlami_lower_bound:TailBoundNegativeInf', ...
                '%s', 'Unexpected: tail bound interval has negative inf. Check interval operations.');
        end
        ddlam_i_upper_bound = R_tail_ub+E_approx;
    end

catch ME
    % MATLAB(MEXCEP)-safe warning with identifier + formatted message
    warning(ME.identifier, '%s', ME.message);
end


fprintf('--- Calculation complete for ddot(lambda_%d) ---\n', i);

fprintf('  - mid( widehat(ddot(lambda_{i,N})) ): %f\n', I_mid(ddlam_i_N_h));
fprintf('  - inf( rigorous lower bound ): %f\n\n', I_inf(ddlam_i_lower_bound));
fprintf('  - sup( tail bound R_N )            : %e\n', I_sup(R_tail_ub));
fprintf('  - sup( rigorous upper bound ): %f\n', I_sup(ddlam_i_upper_bound));




end
