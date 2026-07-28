function mesh = upper_make_reference_mesh(n)
%UPPER_MAKE_REFERENCE_MESH Uniform P1 mesh of conv{(0,0),(1,0),(0,1)}.
%   The integer lattice coordinates are retained so that the certificate
%   evaluator can assemble every reference integral using exact integers
%   and outward-rounded interval arithmetic.

validateattributes(n, {'numeric'}, {'scalar','integer','>=',2});

node_id = zeros(n + 1, n + 1);
grid_ij = zeros((n + 1) * (n + 2) / 2, 2);
node = 0;
for i = 0:n
    for j = 0:(n - i)
        node = node + 1;
        node_id(i + 1, j + 1) = node;
        grid_ij(node, :) = [i, j];
    end
end

triangles = zeros(n * n, 3);
elem = 0;
for i = 0:(n - 1)
    for j = 0:(n - 1 - i)
        elem = elem + 1;
        triangles(elem, :) = [ ...
            node_id(i + 1, j + 1), ...
            node_id(i + 2, j + 1), ...
            node_id(i + 1, j + 2)];
        if i + j <= n - 2
            elem = elem + 1;
            triangles(elem, :) = [ ...
                node_id(i + 2, j + 1), ...
                node_id(i + 2, j + 2), ...
                node_id(i + 1, j + 2)];
        end
    end
end
triangles = triangles(1:elem, :);

boundary = (grid_ij(:, 1) == 0) | (grid_ij(:, 2) == 0) | ...
    (sum(grid_ij, 2) == n);

mesh.n = n;
mesh.grid_ij = grid_ij;
mesh.nodes = grid_ij / n;
mesh.triangles = triangles;
mesh.boundary = boundary;
mesh.interior = find(~boundary);
end
