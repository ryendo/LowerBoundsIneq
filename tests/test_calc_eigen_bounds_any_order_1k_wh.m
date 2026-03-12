



tri_intval = [0, 0, 1, 0, 0.5, sqrt(3)/2];


% Parameters
neig = 3;
N_LG = 16;
N_rho = 64;
LagrangeOrder = 2;


[eig_bounds, LA_eigf, A_grad, A_L2, A_xx, A_xy, A_yy] = calc_eigen_bounds_any_order_1k_wh(neig,tri_intval,N_LG,N_rho,LagrangeOrder);

eig_bounds