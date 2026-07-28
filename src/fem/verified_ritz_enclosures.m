function ritz = verified_ritz_enclosures(U,K,M,neig)
%VERIFIED_RITZ_ENCLOSURES  Rigorous Ritz values on a fixed trial subspace.
%
% Let E_h=span(U(:,1:NEIG)).  The generalized eigenvalues of
%
%   (U'*K*U)c = theta (U'*M*U)c
%
% are the Ritz values on E_h.  Their verified upper endpoints are
% conforming min--max upper bounds for the corresponding exact Dirichlet
% eigenvalues.  In interval mode the compressed problem is solved by
% INTLAB/VEIG; raw floating-point EIGS values are never used as bounds.

global INTERVAL_MODE
if isempty(INTERVAL_MODE)
    INTERVAL_MODE = 0;
end
if nargin < 4 || isempty(neig)
    neig = size(U,2);
end
if neig < 1 || neig ~= floor(neig) || size(U,2) < neig
    error('verified_ritz_enclosures:BadDimension', ...
        'NEIG must be a positive integer not exceeding size(U,2).');
end

V = U(:,1:neig);
KA = V'*K*V;
MB = V'*M*V;
KA = I_hull(KA,KA');
MB = I_hull(MB,MB');

if INTERVAL_MODE
    if ~isspd(MB)
        error('verified_ritz_enclosures:MassNotVerifiedPositiveDefinite', ...
            ['The compressed mass matrix is not verified positive ', ...
             'definite; the supplied trial vectors may be dependent.']);
    end
    [ritz,indices] = veig(KA,MB,1:neig);
    indices = indices(:).';
    if ~isequal(indices,1:neig)
        error('verified_ritz_enclosures:BadVerifiedIndices', ...
            ['VEIG did not return precisely the requested ordered ', ...
             'Ritz indices 1:NEIG.']);
    end
else
    [~,mass_chol_flag] = chol(MB);
    if mass_chol_flag ~= 0
        error('verified_ritz_enclosures:MassNotPositiveDefinite', ...
            ['The compressed mass matrix is not positive definite; ', ...
             'the supplied trial vectors may be dependent.']);
    end
    ritz = eig(KA,MB);
    ritz = sort(real(ritz),'ascend');
end
ritz = ritz(:);
if length(ritz) ~= neig
    error('verified_ritz_enclosures:BadVerifiedShape', ...
        'The compressed generalized eigenproblem returned the wrong size.');
end

for k = 1:neig
    if ~(I_inf(ritz(k)) > 0) || ~isfinite(I_sup(ritz(k)))
        error('verified_ritz_enclosures:InvalidRitzValue', ...
            'Could not verify a finite positive Ritz value %d.',k);
    end
end
end
