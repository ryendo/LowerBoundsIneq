function test_verified_ritz_enclosures(requested_interval_mode)
%TEST_VERIFIED_RITZ_ENCLOSURES  Compressed Ritz/index/cluster regression.

this_file = mfilename('fullpath');
repo_root = fileparts(fileparts(this_file));
addpath(fullfile(repo_root,'src','interval'),'-begin');
addpath(fullfile(repo_root,'src','fem'),'-begin');
addpath(fullfile(repo_root,'src','lib','veigs'),'-begin');

global INTERVAL_MODE
if nargin < 1
    requested_interval_mode = false;
end
INTERVAL_MODE = logical(requested_interval_mode);

% Deliberately descending diagonal order: the helper must return spectral
% indices [1,2], not the storage order inherited from EIG.
K2 = I_intval(diag([4,1]));
M2 = I_eye(2,2);
U2 = I_eye(2,2);
r2 = verified_ritz_enclosures(U2,K2,M2,2);
assert(isequal(size(r2),[2,1]));
if requested_interval_mode
    assert(isa(r2,'intval'));
else
    assert(~isa(r2,'intval'));
end
assert(abs(I_mid(r2(1))-1) <= 1e-12);
assert(abs(I_mid(r2(2))-4) <= 1e-12);
if requested_interval_mode
    assert(I_inf(r2(1)) <= 1 && 1 <= I_sup(r2(1)));
    assert(I_inf(r2(2)) <= 4 && 4 <= I_sup(r2(2)));
end
assert(I_sup(r2(1)) <= I_sup(r2(2)));

% Nonzero coefficient widths exercise actual interval-matrix
% certification rather than merely roundoff enclosure of a point problem.
Kwide = I_infsup(diag([0.9,3.9]),diag([1.1,4.1]));
rwide = verified_ritz_enclosures(U2,Kwide,M2,2);
assert(isequal(size(rwide),[2,1]));
if requested_interval_mode
    assert(I_inf(rwide(1)) <= 0.9 && 1.1 <= I_sup(rwide(1)));
    assert(I_inf(rwide(2)) <= 3.9 && 4.1 <= I_sup(rwide(2)));
else
    assert(all(abs(rwide-[1;4]) <= 1e-12));
end

% A multiple compressed eigenvalue must retain both spectral indices.
Krep = I_intval(2*eye(2));
rrep = verified_ritz_enclosures(U2,Krep,M2,2);
assert(isequal(size(rrep),[2,1]));
assert(all(abs(I_mid(rrep)-2) <= 1e-12));
if requested_interval_mode
    assert(all(I_inf(rrep) <= 2) && all(2 <= I_sup(rrep)));
end

r1 = verified_ritz_enclosures(U2(:,1),K2,M2,1);
assert(isequal(size(r1),[1,1]));
assert(abs(I_mid(r1)-4) <= 1e-12);
if requested_interval_mode
    assert(I_inf(r1) <= 4 && 4 <= I_sup(r1));
end

% A dependent trial basis must be rejected before invoking the verified
% generalized eigenvalue solver, whose theory requires a positive
% definite compressed mass matrix.
caught_dependent_basis = false;
try
    verified_ritz_enclosures([U2(:,1),U2(:,1)],K2,M2,2);
catch ME
    caught_dependent_basis = startsWith( ...
        ME.identifier,'verified_ritz_enclosures:MassNot');
end
assert(caught_dependent_basis);

% Cross terms make the maximum Rayleigh quotient on a two-dimensional
% cluster strictly larger than both tested column quotients.
K3 = I_intval(diag([1,4,9]));
M3 = I_eye(3,3);
U3 = I_intval([1,1,0; 1,-0.2,0; 0,0,1]);
individual = I_intval(zeros(3,1));
for k = 1:3
    individual(k) = ...
        (U3(:,k)'*K3*U3(:,k))/(U3(:,k)'*M3*U3(:,k));
end
cluster_ritz = verified_ritz_enclosures( ...
    U3(:,1:2),K3,M3,2);
assert(abs(I_mid(cluster_ritz(2))-4) <= 1e-12);
if requested_interval_mode
    assert(I_inf(cluster_ritz(2)) <= 4 ...
        && 4 <= I_sup(cluster_ritz(2)));
end
assert(I_inf(cluster_ritz(2)) > max(I_sup(individual(1:2))));

lambda = I_intval([1;4;9]);
clusters = {[1,2],[3]};
[~,~,~,info] = calc_grad_error_bounds( ...
    lambda,individual,U3,K3,M3,clusters);
assert(info.ok(1));
assert(abs(I_mid(info.lambda_Nkh(1))-4) <= 1e-12);
if requested_interval_mode
    % lambda_Nkh is deliberately stored as the point interval at the
    % verified upper endpoint, rather than as a two-sided Ritz enclosure.
    assert(I_inf(info.lambda_Nkh(1)) >= 4);
end
assert(I_inf(info.lambda_Nkh(1)) ...
    > max(I_sup(individual(1:2))));

fprintf('test_verified_ritz_enclosures: PASS (INTERVAL_MODE=%d)\n', ...
    INTERVAL_MODE);
end
