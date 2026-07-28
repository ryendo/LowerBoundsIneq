function mesh = apply_exact_boundary_point_setting(mesh,tri_vertices)
%APPLY_EXACT_BOUNDARY_POINT_SETTING  Certified triangle-boundary snapping.
%
% Every topological boundary node is either snapped to a triangle vertex or
% written constructively as
%
%        p_j + t (p_k-p_j),             0 <= t <= 1,
%
% on one of the three domain sides.  Classification is performed using the
% midpoint Gmsh geometry, whereas the construction uses TRI_VERTICES itself;
% hence interval vertices produce an interval enclosure on the corresponding
% moving side.  Any unclassifiable node is an error: no unchanged boundary
% coordinate is allowed to escape certification.
%
% A fail-closed postcheck proves:
%   (i) every topological boundary node has the construction above;
%  (ii) the endpoints of every boundary mesh edge share a domain side;
% (iii) every interval element signed area excludes zero and has the same
%       orientation as TRI_VERTICES.

validate_mesh_input(mesh,tri_vertices);

nodes = I_intval(mesh.nodes);
bedges = double(mesh.boundary_edges);
elements = double(mesh.elements);
vertices = I_intval(tri_vertices);
mid_vertices = I_mid(vertices);

domain_area2 = signed_area_twice(vertices(1,:),vertices(2,:),vertices(3,:));
domain_sign = certified_nonzero_sign(domain_area2, ...
    'apply_exact_boundary_point_setting:DegenerateTriangle', ...
    'The interval signed area of TRI_VERTICES contains zero.');

% Gmsh writes the midpoint geometry in double precision.  The original
% absolute tolerance is retained, but failure now aborts certification.
vertex_tol = 1e-6;
edge_tol = 1e-6;

number_of_nodes = size(nodes,1);
constructed = false(number_of_nodes,1);
construction_kind = zeros(number_of_nodes,1); % 1: vertex, 2: side
construction_id = zeros(number_of_nodes,1);
side_membership = false(number_of_nodes,3);
side_parameter = NaN(number_of_nodes,3);

boundary_nodes = unique(bedges(:))';
side_start = [1,2,3];
side_end = [2,3,1];

for idx = boundary_nodes
    x = I_mid(nodes(idx,:));

    vertex_distances = [ ...
        norm(x-mid_vertices(1,:)), ...
        norm(x-mid_vertices(2,:)), ...
        norm(x-mid_vertices(3,:))];
    [nearest_vertex_distance,vertex_id] = min(vertex_distances);

    if nearest_vertex_distance <= vertex_tol
        nodes(idx,:) = vertices(vertex_id,:);
        constructed(idx) = true;
        construction_kind(idx) = 1;
        construction_id(idx) = vertex_id;
        [side_membership(idx,:),side_parameter(idx,:)] = ...
            vertex_side_data(vertex_id);
        continue
    end

    distances = zeros(1,3);
    parameters = zeros(1,3);
    for side_id = 1:3
        [parameters(side_id),distances(side_id)] = ...
            point_to_segment_projection_data( ...
                x,mid_vertices(side_start(side_id),:), ...
                mid_vertices(side_end(side_id),:));
    end
    [nearest_side_distance,side_id] = min(distances);
    if ~(isfinite(nearest_side_distance) && nearest_side_distance <= edge_tol)
        error('apply_exact_boundary_point_setting:BoundaryNodeClassification', ...
            ['Boundary node %d is %.17g from its nearest triangle side; ' ...
             'the admissible distance is %.17g.'], ...
            idx,nearest_side_distance,edge_tol);
    end

    t = parameters(side_id);
    if ~(isfinite(t) && t >= 0 && t <= 1)
        error('apply_exact_boundary_point_setting:InvalidProjectionParameter', ...
            'Boundary node %d produced the invalid side parameter %.17g.', ...
            idx,t);
    end
    a = vertices(side_start(side_id),:);
    b = vertices(side_end(side_id),:);
    nodes(idx,:) = a+t*(b-a);
    constructed(idx) = true;
    construction_kind(idx) = 2;
    construction_id(idx) = side_id;
    side_membership(idx,side_id) = true;
    side_parameter(idx,side_id) = t;
end

postcheck_boundary_construction(nodes,vertices,boundary_nodes, ...
    constructed,construction_kind,construction_id,side_membership, ...
    side_parameter);
postcheck_boundary_edges(bedges,elements,side_membership,side_parameter);
postcheck_element_orientations(nodes,elements,domain_sign);

mesh.nodes = nodes;
end


function validate_mesh_input(mesh,tri_vertices)
required_fields = {'nodes','boundary_edges','elements'};
if ~isstruct(mesh)
    error('apply_exact_boundary_point_setting:InvalidMesh', ...
        'MESH must be a structure.');
end
for k = 1:numel(required_fields)
    if ~isfield(mesh,required_fields{k})
        error('apply_exact_boundary_point_setting:InvalidMesh', ...
            'MESH is missing the field "%s".',required_fields{k});
    end
end
if size(mesh.nodes,2) ~= 2 || isempty(mesh.nodes)
    error('apply_exact_boundary_point_setting:InvalidMesh', ...
        'MESH.nodes must be a nonempty N-by-2 array.');
end
if ~isequal(size(tri_vertices),[3,2])
    error('apply_exact_boundary_point_setting:InvalidTriangle', ...
        'TRI_VERTICES must be a 3-by-2 array.');
end
if size(mesh.boundary_edges,2) ~= 2 || isempty(mesh.boundary_edges)
    error('apply_exact_boundary_point_setting:InvalidMesh', ...
        'MESH.boundary_edges must be a nonempty M-by-2 array.');
end
if size(mesh.elements,2) ~= 3 || isempty(mesh.elements)
    error('apply_exact_boundary_point_setting:InvalidMesh', ...
        'MESH.elements must be a nonempty K-by-3 array.');
end
interval_nodes = I_intval(mesh.nodes(:));
interval_vertices = I_intval(tri_vertices(:));
if any(~isfinite(I_inf(interval_nodes))) || ...
        any(~isfinite(I_sup(interval_nodes))) || ...
        any(~isfinite(I_inf(interval_vertices))) || ...
        any(~isfinite(I_sup(interval_vertices)))
    error('apply_exact_boundary_point_setting:NonfiniteCoordinate', ...
        'All node and triangle-vertex bounds must be finite.');
end

number_of_nodes = size(mesh.nodes,1);
index_arrays = {mesh.boundary_edges,mesh.elements};
for k = 1:numel(index_arrays)
    indices = double(index_arrays{k}(:));
    if any(~isfinite(indices)) || any(indices ~= floor(indices)) || ...
            any(indices < 1) || any(indices > number_of_nodes)
        error('apply_exact_boundary_point_setting:InvalidConnectivity', ...
            'All connectivity entries must be valid integer node indices.');
    end
end
if any(mesh.boundary_edges(:,1) == mesh.boundary_edges(:,2))
    error('apply_exact_boundary_point_setting:InvalidBoundaryEdge', ...
        'A boundary edge has identical endpoint indices.');
end
if any(mesh.elements(:,1) == mesh.elements(:,2) | ...
       mesh.elements(:,2) == mesh.elements(:,3) | ...
       mesh.elements(:,3) == mesh.elements(:,1))
    error('apply_exact_boundary_point_setting:InvalidElement', ...
        'A triangle element repeats a node index.');
end
end


function [membership,parameters] = vertex_side_data(vertex_id)
membership = false(1,3);
parameters = NaN(1,3);
switch vertex_id
    case 1
        membership([1,3]) = true;
        parameters([1,3]) = [0,1];
    case 2
        membership([1,2]) = true;
        parameters([1,2]) = [1,0];
    case 3
        membership([2,3]) = true;
        parameters([2,3]) = [1,0];
    otherwise
        error('apply_exact_boundary_point_setting:InternalVertexIndex', ...
            'Unexpected triangle vertex index %d.',vertex_id);
end
end


function [t,distance] = point_to_segment_projection_data(point,a,b)
direction = b-a;
length_squared = dot(direction,direction);
if ~(isfinite(length_squared) && length_squared > 0)
    error('apply_exact_boundary_point_setting:DegenerateMidpointSide', ...
        'A midpoint triangle side has zero or nonfinite length.');
end
t = dot(point-a,direction)/length_squared;
t = max(0,min(1,t));
projection = a+t*direction;
distance = norm(point-projection);
end


function postcheck_boundary_construction(nodes,vertices,boundary_nodes, ...
        constructed,construction_kind,construction_id,side_membership, ...
        side_parameter)
side_start = [1,2,3];
side_end = [2,3,1];

for idx = boundary_nodes
    if ~constructed(idx) || ~any(side_membership(idx,:))
        error('apply_exact_boundary_point_setting:UnconstructedBoundaryNode', ...
            'Boundary node %d has no certified triangle-side construction.',idx);
    end

    member_parameters = side_parameter(idx,side_membership(idx,:));
    if any(~isfinite(member_parameters)) || ...
            any(member_parameters < 0) || any(member_parameters > 1)
        error('apply_exact_boundary_point_setting:InvalidBoundaryCertificate', ...
            'Boundary node %d has invalid side-membership parameters.',idx);
    end

    if construction_kind(idx) == 1
        expected = vertices(construction_id(idx),:);
    elseif construction_kind(idx) == 2
        side_id = construction_id(idx);
        t = side_parameter(idx,side_id);
        a = vertices(side_start(side_id),:);
        b = vertices(side_end(side_id),:);
        expected = a+t*(b-a);
    else
        error('apply_exact_boundary_point_setting:InvalidBoundaryCertificate', ...
            'Boundary node %d has no recognized construction type.',idx);
    end

    if ~same_interval_vector(nodes(idx,:),expected)
        error('apply_exact_boundary_point_setting:BoundaryConstructionMismatch', ...
            'Boundary node %d does not equal its certified construction.',idx);
    end
end
end


function postcheck_boundary_edges(bedges,elements,side_membership,side_parameter)
side_used = false(1,3);
for k = 1:size(bedges,1)
    i = bedges(k,1);
    j = bedges(k,2);
    common_sides = find(side_membership(i,:) & side_membership(j,:));
    if isempty(common_sides)
        error('apply_exact_boundary_point_setting:BoundaryEdgeDomainSide', ...
            ['Boundary edge (%d,%d) has endpoints on different triangle ' ...
             'sides.'],i,j);
    end
    parameter_difference = abs( ...
        side_parameter(i,common_sides)-side_parameter(j,common_sides));
    nondegenerate_sides = common_sides(parameter_difference > 0);
    if isempty(nondegenerate_sides)
        error('apply_exact_boundary_point_setting:CollapsedBoundaryEdge', ...
            'Boundary edge (%d,%d) collapses under boundary projection.',i,j);
    end
    side_used(nondegenerate_sides) = true;
end
if ~all(side_used)
    error('apply_exact_boundary_point_setting:IncompleteDomainBoundary', ...
        'The mesh boundary does not constructively cover all triangle sides.');
end

declared_edges = sort(double(bedges),2);
unique_declared_edges = unique(declared_edges,'rows');
if size(unique_declared_edges,1) ~= size(declared_edges,1)
    error('apply_exact_boundary_point_setting:DuplicateBoundaryEdge', ...
        'MESH.boundary_edges contains duplicate undirected edges.');
end

all_element_edges = sort([ ...
    elements(:,[1,2]); ...
    elements(:,[2,3]); ...
    elements(:,[3,1])],2);
[unique_element_edges,~,edge_map] = unique(all_element_edges,'rows');
edge_counts = accumarray(edge_map,1);
if any(edge_counts > 2)
    error('apply_exact_boundary_point_setting:NonmanifoldMesh', ...
        'An element edge has incidence greater than two.');
end
topological_boundary = unique_element_edges(edge_counts == 1,:);
if ~isequal(sortrows(unique_declared_edges),sortrows(topological_boundary))
    error('apply_exact_boundary_point_setting:BoundaryTopologyMismatch', ...
        ['MESH.boundary_edges does not equal the boundary induced by ' ...
         'MESH.elements.']);
end
end


function postcheck_element_orientations(nodes,elements,domain_sign)
for k = 1:size(elements,1)
    ids = elements(k,:);
    area2 = signed_area_twice( ...
        nodes(ids(1),:),nodes(ids(2),:),nodes(ids(3),:));
    element_sign = certified_nonzero_sign(area2, ...
        'apply_exact_boundary_point_setting:ElementAreaContainsZero', ...
        sprintf(['The interval signed area of element %d contains zero ' ...
                 'or is nonfinite.'],k));
    if element_sign ~= domain_sign
        error('apply_exact_boundary_point_setting:InconsistentElementOrientation', ...
            ['Element %d has orientation opposite to TRI_VERTICES; all ' ...
             'element orientations must agree.'],k);
    end
end
end


function area2 = signed_area_twice(a,b,c)
ab = b-a;
ac = c-a;
area2 = ab(1)*ac(2)-ab(2)*ac(1);
end


function sign_value = certified_nonzero_sign(value,error_id,error_message)
lower = I_inf(value);
upper = I_sup(value);
if ~(isscalar(lower) && isscalar(upper) && ...
        isfinite(lower) && isfinite(upper))
    error(error_id,'%s',error_message);
end
if lower > 0
    sign_value = 1;
elseif upper < 0
    sign_value = -1;
else
    error(error_id,'%s',error_message);
end
end


function tf = same_interval_vector(a,b)
tf = isequal(I_inf(a),I_inf(b)) && isequal(I_sup(a),I_sup(b));
end
