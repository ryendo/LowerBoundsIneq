


neig = 3;
lagrange_order = 1;

mesh = make_mesh_by_gmsh(0.5, sqrt(3)/2, 0.1);
vert = I_intval(mesh.nodes);
edge = mesh.edges;
tri  = mesh.elements;
bd   = mesh.boundary_edges;


[eig_value, eig_func_no_bdry, eig_func_with_bdry, A_grad, A_L2, A_xx, A_xy, A_yy, bd_dof_idx] = laplace_eig_lagrange_detailed(lagrange_order, vert, edge, tri, bd, neig);
