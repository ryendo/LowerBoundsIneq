function test_omega_up_all_functionals_unified()
%TEST_OMEGA_UP_ALL_FUNCTIONALS_UNIFIED  Geometry reuse and disk-skip test.

this_file = mfilename('fullpath');
repo_root = fileparts(fileparts(this_file));
addpath(fullfile(repo_root,'src','interval'),'-begin');
addpath(fullfile(repo_root,'src','algorithms'),'-begin');

global INTERVAL_MODE
INTERVAL_MODE = 0;

x = 0.53;
y = 0.84;
dlambda = [-7,-3]; %#ok<NASGU>

% In x direction dlambda is irrelevant.  Changing only the Hessian upper
% endpoint must alter JUP but leave both lower-functional bounds unchanged.
b1 = omega_up_all_functional_curvatures('x',x,y, ...
    -5,20,30);
b2 = omega_up_all_functional_curvatures('x',x,y, ...
    100,20,40);
assert(abs(b1.J1.bound-b2.J1.bound) < 1e-12);
assert(abs(b1.J2.bound-b2.J2.bound) < 1e-12);
assert(abs((b2.JUP.bound-b1.JUP.bound)-0.5*y*10) < 1e-12);

% In y direction the same first-derivative enclosure is shared by all
% functionals.  A scalar shift of dlambda translates every bound equally.
by1 = omega_up_all_functional_curvatures('y',x,y, ...
    -5,20,30);
by2 = omega_up_all_functional_curvatures('y',x,y, ...
    -2,20,30);
assert(abs((by2.J1.bound-by1.J1.bound)-3) < 1e-12);
assert(abs((by2.J2.bound-by1.J2.bound)-3) < 1e-12);
assert(abs((by2.JUP.bound-by1.JUP.bound)-3) < 1e-12);

% An infinite signed endpoint is not a certificate, even when its formal
% comparison with zero has the desired sign.
bad_lower = omega_up_all_functional_curvatures( ...
    'x',x,y,-5,Inf,30);
bad_upper = omega_up_all_functional_curvatures( ...
    'x',x,y,-5,20,-Inf);
assert(~bad_lower.J1.sign_verified);
assert(~bad_lower.J2.sign_verified);
assert(~bad_upper.JUP.sign_verified);

% A rectangle whose lower-left corner lies outside the unit disk is
% coordinatewise outside.  It must be skipped before any FEM/RT call.
cell_def = struct( ...
    'id',1,'kind','rectangle','ix',1,'iy',1, ...
    'x_lo',0.74,'x_hi',0.741,'y_lo',0.75,'y_hi',0.751);
mesh_params = struct( ...
    'N_LG',4,'N_rho',12,'bound_order',2, ...
    'trial_order',3,'N_trial',4,'RT_order',3);
record = verify_omega_up_all_functionals_cell(cell_def,mesh_params);
assert(strcmp(record.status,'skipped_outside_disk'));
assert(~record.applicable && ~record.computed);
assert(record.estimator_calls == 0);

% Functional scope is independent of geometric applicability.  A cell with
% no requested functional is a neutral out-of-scope record and must also be
% skipped before any FEM/RT call.
cell_def = struct( ...
    'id',2,'kind','axis','ix',0,'iy',1, ...
    'x_lo',0.5,'x_hi',0.5,'y_lo',0.84,'y_hi',0.841, ...
    'J1_in_scope',false,'J2_in_scope',false,'JUP_in_scope',false);
record = verify_omega_up_all_functionals_cell(cell_def,mesh_params);
assert(strcmp(record.status,'skipped_out_of_scope'));
assert(~record.applicable && record.domain_applicable && ~record.computed);
assert(record.estimator_calls == 0);
assert(record.all_scoped_signs_ok);
assert(~record.J1_certified && ~record.J2_certified ...
    && ~record.JUP_certified);

fprintf('test_omega_up_all_functionals_unified: PASS\n');
end
