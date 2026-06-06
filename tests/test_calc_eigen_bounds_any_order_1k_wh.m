



tri_intval = [0, 0, 1, 0, 0.5, sqrt(3)/2];


% Parameters
neig = 1;
N_LG = 8;
N_rho = 32;
LagrangeOrder = 3;
fem_ord_LG = 4;


[eig_bounds, LA_eigf, A_grad, A_L2, A_xx, A_xy, A_yy] = calc_eigen_bounds_any_order_1k_wh(neig,tri_intval,N_LG,N_rho,LagrangeOrder);

eig_bounds