function invs = upper_trial_invariants(mesh, coeff, mode)
%UPPER_TRIAL_INVARIANTS Certified scalar integrals for a fixed P1 trial.
%   In interval mode the binary64 coefficients are point intervals.  The
%   reference geometry uses integer lattice data and every division is
%   performed by INTLAB, so no rounded sparse matrix is trusted.

if nargin < 3
    mode = 'double';
end
is_interval = strcmpi(mode, 'interval');
if is_interval
    zero = intval(0);
    area_scale = intval(1) / (2 * mesh.n^2);
else
    zero = 0;
    area_scale = 1 / (2 * mesh.n^2);
end

Exx = zero;
Exy = zero;
Eyy = zero;
M0 = zero;
for t = 1:size(mesh.triangles,1)
    ids = mesh.triangles(t,:);
    ij = mesh.grid_ij(ids,:);
    det_int = (ij(2,1)-ij(1,1)) * (ij(3,2)-ij(1,2)) - ...
              (ij(3,1)-ij(1,1)) * (ij(2,2)-ij(1,2));
    area = area_scale * abs(det_int);
    bx = [ij(2,2)-ij(3,2); ij(3,2)-ij(1,2); ij(1,2)-ij(2,2)];
    by = [ij(3,1)-ij(2,1); ij(1,1)-ij(3,1); ij(2,1)-ij(1,1)];
    if is_interval
        q = intval(coeff(ids));
        gx = intval(mesh.n) / intval(det_int) * (intval(bx') * q);
        gy = intval(mesh.n) / intval(det_int) * (intval(by') * q);
    else
        q = coeff(ids);
        gx = mesh.n / det_int * (bx' * q);
        gy = mesh.n / det_int * (by' * q);
    end
    Exx = Exx + area * gx^2;
    Exy = Exy + area * gx*gy;
    Eyy = Eyy + area * gy^2;
    M0 = M0 + area/6 * ( ...
        q(1)^2 + q(2)^2 + q(3)^2 + ...
        q(1)*q(2) + q(1)*q(3) + q(2)*q(3));
end

invs.Exx = Exx;
invs.Exy = Exy;
invs.Eyy = Eyy;
invs.M0 = M0;
end
