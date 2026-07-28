function model = upper_make_high_order_model(trial_order, N_trial, mode)
%UPPER_MAKE_HIGH_ORDER_MODEL Build one affine-transportable Pp trial space.
%   A single Gmsh mesh on T(1/2,1/2) is assembled with the existing VFEM
%   Lagrange code. Its four matrices are pulled back to the standard
%   reference triangle, so all atlas cells reuse one topology and Gmsh is
%   not called per cell.

if nargin < 3
    mode = 'interval';
end
validateattributes(trial_order, {'numeric'}, {'scalar','integer','>=',1});
validateattributes(N_trial, {'numeric'}, {'scalar','integer','>=',2});

x_ref = 0.5;
y_ref = 0.5;
mesh = make_mesh_by_gmsh(I_intval(x_ref), I_intval(y_ref), 1/N_trial, ...
    [], false, sprintf('upper-p%d-n%d', trial_order, N_trial));
[~,~,~,~,M_phys,Kxx_phys,Kxy_phys,Kyy_phys] = ...
    laplace_eig_lagrange_detailed( ...
        trial_order, mesh.nodes, mesh.edges, mesh.elements, ...
        mesh.boundary_edges, 1);

% If B=[1 x_ref;0 y_ref], then
% C_ref=B'*C_phys*B/det(B), M_ref=M_phys/det(B).
Kxx_ref = Kxx_phys/y_ref;
Kxy_ref = (x_ref*Kxx_phys+y_ref*Kxy_phys)/y_ref;
Kyy_ref = (x_ref^2*Kxx_phys+2*x_ref*y_ref*Kxy_phys ...
    +y_ref^2*Kyy_phys)/y_ref;
M_ref = M_phys/y_ref;

model.backend = 'lagrange_high_order';
model.mode = lower(mode);
model.trial_order = trial_order;
model.N_trial = N_trial;
model.reference_x = x_ref;
model.reference_y = y_ref;
model.Kxx = Kxx_ref;
model.Kxy = Kxy_ref;
model.Kyy = Kyy_ref;
model.M = M_ref;
model.ndof = size(M_ref,1);
end
