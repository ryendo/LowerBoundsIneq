function mesh = affine_map_triangle_mesh(mesh_ref, a0_mid, b0_mid, a, b)
% Map mesh on triangle (0,0)-(1,0)-(a0_mid,b0_mid) to (0,0)-(1,0)-(a,b)

    mesh = mesh_ref;

    % shear + scale parameters
    alpha = (a - a0_mid) / b0_mid;
    beta  = b / b0_mid;

    xy = mesh_ref.nodes;  % typically intval after apply_exact_boundary_point_setting

    % (x,y) -> (x + alpha*y, beta*y)
    mesh.nodes = [ xy(:,1) + alpha .* xy(:,2),  beta .* xy(:,2) ];

    % update domain
    mesh.domain = I_intval([0 0; 1 0; a b]);
end