function [coeff, lambda_ritz] = upper_midpoint_trial(mesh, matrices, x, y)
%UPPER_MIDPOINT_TRIAL First conforming P1 Ritz vector at T(x,y).
%   This solve is only a coefficient generator.  Certification never uses
%   lambda_ritz: the returned binary64 coefficients are subsequently held
%   fixed and their Rayleigh quotient is re-evaluated with INTLAB.

if ~(isfinite(x) && isfinite(y) && y > 0)
    error('A nondegenerate finite triangle is required.');
end

ii = mesh.interior;
K0 = (y^2 + x^2) * matrices.Kxx(ii,ii) ...
    - 2*x * matrices.Kxy(ii,ii) + matrices.Kyy(ii,ii);
M0 = y^2 * matrices.M(ii,ii);

try
    if exist('OCTAVE_VERSION', 'builtin')
        opts.tol = 1e-11;
        opts.maxit = 2000;
        [v, d] = eigs(K0, M0, 1, 'sm', opts);
    else
        opts.tol = 1e-11;
        opts.maxit = 2000;
        opts.isreal = true;
        opts.issym = true;
        [v, d] = eigs(K0, M0, 1, 'smallestreal', opts);
    end
    lambda_ritz = real(d(1,1));
catch
    [V, D] = eig(full(K0), full(M0));
    values = real(diag(D));
    values(values <= 0 | ~isfinite(values)) = Inf;
    [lambda_ritz, idx] = min(values);
    v = V(:, idx);
end

v = real(v);
v = v / max(abs(v));
coeff = zeros(size(mesh.grid_ij,1), 1);
coeff(ii) = v;
end
