function test_apply_exact_boundary_point_setting(use_current_mode)
%TEST_APPLY_EXACT_BOUNDARY_POINT_SETTING  Fail-closed boundary regression.

this_file = mfilename('fullpath');
repo_root = fileparts(fileparts(this_file));
addpath(fullfile(repo_root,'src','interval'),'-begin');
addpath(fullfile(repo_root,'src','mesh'),'-begin');

global INTERVAL_MODE
if nargin < 1 || ~logical(use_current_mode)
    INTERVAL_MODE = 0;
end

vertices = I_intval([0,0;1,0;0.3,0.8]);
if logical(INTERVAL_MODE)
    % The projected nodes must follow an interval-valued moving side, not
    % merely the midpoint side used for Gmsh classification.
    vertices(3,1) = infsup(0.29999,0.30001);
    vertices(3,2) = infsup(0.79999,0.80001);
end
mesh = valid_noisy_mesh();
certified_mesh = apply_exact_boundary_point_setting(mesh,vertices);

% Vertex snapping is exact, and every projected side point satisfies its
% defining affine line in the point (noninterval) regression.
assert(same_bounds(certified_mesh.nodes(1,:),vertices(1,:)));
assert(same_bounds(certified_mesh.nodes(3,:),vertices(2,:)));
assert(same_bounds(certified_mesh.nodes(5,:),vertices(3,:)));
if ~logical(INTERVAL_MODE)
    assert(abs(certified_mesh.nodes(2,2)) < 1e-15);
    assert(abs(cross2(certified_mesh.nodes(4,:)-vertices(2,:), ...
        vertices(3,:)-vertices(2,:))) < 1e-15);
    assert(abs(cross2(certified_mesh.nodes(6,:)-vertices(3,:), ...
        vertices(1,:)-vertices(3,:))) < 1e-15);
end

bad = mesh;
bad.nodes(2,:) = [0.5,2e-6];
assert_throws(@() apply_exact_boundary_point_setting(bad,vertices), ...
    'apply_exact_boundary_point_setting:BoundaryNodeClassification');

bad = mesh;
bad.boundary_edges(1,:) = [2,4];
assert_throws(@() apply_exact_boundary_point_setting(bad,vertices), ...
    'apply_exact_boundary_point_setting:BoundaryEdgeDomainSide');

bad = mesh;
bad.boundary_edges(1,:) = [];
assert_throws(@() apply_exact_boundary_point_setting(bad,vertices), ...
    'apply_exact_boundary_point_setting:BoundaryTopologyMismatch');

bad = mesh;
bad.elements(1,[1,2]) = bad.elements(1,[2,1]);
assert_throws(@() apply_exact_boundary_point_setting(bad,vertices), ...
    'apply_exact_boundary_point_setting:InconsistentElementOrientation');

if logical(INTERVAL_MODE)
    bad = mesh;
    bad.nodes = I_intval(bad.nodes);
    bad.nodes(7,1) = infsup(-1,2);
    assert_throws(@() apply_exact_boundary_point_setting(bad,vertices), ...
        'apply_exact_boundary_point_setting:ElementAreaContainsZero');
end

fprintf('test_apply_exact_boundary_point_setting: PASS (interval=%d)\n', ...
    logical(INTERVAL_MODE));
end


function mesh = valid_noisy_mesh()
mesh.nodes = [ ...
     2e-8,-3e-8; ...
     0.5,2e-7; ...
     1-2e-8,3e-8; ...
     0.6500001,0.4000001; ...
     0.3-2e-8,0.8+2e-8; ...
     0.1499999,0.3999999; ...
     0.433333333333333,0.266666666666667];
mesh.boundary_edges = [1,2;2,3;3,4;4,5;5,6;6,1];
mesh.elements = [ ...
    1,2,7; ...
    2,3,7; ...
    3,4,7; ...
    4,5,7; ...
    5,6,7; ...
    6,1,7];
end


function value = cross2(a,b)
value = a(1)*b(2)-a(2)*b(1);
end


function tf = same_bounds(a,b)
tf = isequal(I_inf(a),I_inf(b)) && isequal(I_sup(a),I_sup(b));
end


function assert_throws(action,expected_identifier)
try
    action();
catch exception
    assert(strcmp(exception.identifier,expected_identifier), ...
        'Expected "%s", received "%s": %s', ...
        expected_identifier,exception.identifier,exception.message);
    return
end
error('test_apply_exact_boundary_point_setting:MissingError', ...
    'Expected the error "%s".',expected_identifier);
end
