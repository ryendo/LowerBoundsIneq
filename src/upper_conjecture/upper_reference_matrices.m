function matrices = upper_reference_matrices(mesh)
%UPPER_REFERENCE_MATRICES Assemble P1 derivative and mass matrices.
%   For a coefficient vector c,
%     c'*Kxx*c = int u_xi^2,
%     c'*Kxy*c = int u_xi*u_eta,
%     c'*Kyy*c = int u_eta^2,
%     c'*M*c   = int u^2
%   on the reference right triangle.

nt = size(mesh.triangles, 1);
nn = size(mesh.grid_ij, 1);
rows = zeros(9 * nt, 1);
cols = zeros(9 * nt, 1);
vxx = zeros(9 * nt, 1);
vxy = zeros(9 * nt, 1);
vyy = zeros(9 * nt, 1);
vm = zeros(9 * nt, 1);
cursor = 0;

for t = 1:nt
    ids = mesh.triangles(t, :);
    ij = mesh.grid_ij(ids, :);
    det_int = (ij(2,1)-ij(1,1)) * (ij(3,2)-ij(1,2)) - ...
              (ij(3,1)-ij(1,1)) * (ij(2,2)-ij(1,2));
    if det_int == 0
        error('Degenerate reference element %d.', t);
    end
    area = abs(det_int) / (2 * mesh.n^2);
    gx = mesh.n / det_int * [ ...
        ij(2,2)-ij(3,2); ij(3,2)-ij(1,2); ij(1,2)-ij(2,2)];
    gy = mesh.n / det_int * [ ...
        ij(3,1)-ij(2,1); ij(1,1)-ij(3,1); ij(2,1)-ij(1,1)];
    local_xx = area * (gx * gx');
    local_xy = area * 0.5 * (gx * gy' + gy * gx');
    local_yy = area * (gy * gy');
    local_m = area / 12 * (ones(3) + eye(3));
    for a = 1:3
        for b = 1:3
            cursor = cursor + 1;
            rows(cursor) = ids(a);
            cols(cursor) = ids(b);
            vxx(cursor) = local_xx(a,b);
            vxy(cursor) = local_xy(a,b);
            vyy(cursor) = local_yy(a,b);
            vm(cursor) = local_m(a,b);
        end
    end
end

matrices.Kxx = sparse(rows, cols, vxx, nn, nn);
matrices.Kxy = sparse(rows, cols, vxy, nn, nn);
matrices.Kyy = sparse(rows, cols, vyy, nn, nn);
matrices.M = sparse(rows, cols, vm, nn, nn);
end
