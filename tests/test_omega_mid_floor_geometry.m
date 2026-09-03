function test_omega_mid_floor_geometry()
%TEST_OMEGA_MID_FLOOR_GEOMETRY  Regression for the y=0.04 comparison.

global INTERVAL_MODE
if isempty(INTERVAL_MODE) || ~logical(INTERVAL_MODE)
    error('test_omega_mid_floor_geometry:IntervalModeRequired', ...
        'Initialize INTLAB interval mode before this test.');
end

project_root = fileparts(fileparts(mfilename('fullpath')));
input_table = local_read_input_table( ...
    fullfile(project_root,'inputs','cell_def.csv'));

% Historical input cell 1363 crosses y=0.04.  Read the endpoint strings
% exactly as stored in the proof input, without a binary-double round trip.
cell1363 = local_input_cell(input_table,1363);
assert(strcmp(cell1363.x_inf,'0.53515999999999997'));
assert(strcmp(cell1363.x_sup,'0.53905999999999998'));
assert(strcmp(cell1363.theta_inf,'0.073631000000000002'));
assert(strcmp(cell1363.theta_sup,'0.081811999999999996'));
[area1363,perimeter1363,diagnostics1363] = ...
    local_geometry(cell1363);
local_assert_floor_comparison( ...
    area1363,diagnostics1363,cell1363.x_sup);
local_assert_target_in_cell(cell1363,'0.537','0.041');
local_assert_functional_comparison( ...
    area1363,perimeter1363,'0.537','0.041');

% Cell 3148 is the important endpoint case x_sup=1.  Its comparison point
% (1,0.04) lies outside Omega, as does its upper spectral corner, but both
% are legitimate auxiliary triangles and the B_1/B_2 comparisons hold.
cell3148 = local_input_cell(input_table,3148);
assert(strcmp(cell3148.x_inf,'0.99609499999999995'));
assert(strcmp(cell3148.x_sup,'1'));
assert(strcmp(cell3148.theta_inf,'0.032724999999999997'));
assert(strcmp(cell3148.theta_sup,'0.040905999999999998'));
[area3148,perimeter3148,diagnostics3148] = ...
    local_geometry(cell3148);
local_assert_floor_comparison( ...
    area3148,diagnostics3148,cell3148.x_sup);
one = I_intval('1');
q_norm_squared = diagnostics3148.p_comparison.x^2 ...
    +diagnostics3148.p_comparison.y^2;
upper_norm_squared = diagnostics3148.p_upper.x^2 ...
    +diagnostics3148.p_upper.y^2;
assert(I_inf(q_norm_squared) > I_sup(one));
assert(I_inf(upper_norm_squared) > I_sup(one));
local_assert_target_in_cell(cell3148,'0.999','0.0405');
local_assert_functional_comparison( ...
    area3148,perimeter3148,'0.999','0.0405');

% The former clipping rule (x_inf,0.04) is an overestimate.  At the same
% height, a target point farther right has a larger perimeter and hence a
% strictly smaller B_k.  The corrected x_sup comparison remains below it.
y_floor = I_intval('0.04');
x_left = I_intval(cell3148.x_inf);
old_area = y_floor/I_intval('2');
old_perimeter = one+sqrt(x_left^2+y_floor^2) ...
    +sqrt((one-x_left)^2+y_floor^2);
local_assert_target_in_cell(cell3148,'0.999','0.04',false);
for kind = {'J1','J2'}
    old_value = local_B_from_geometry( ...
        kind{1},old_area,old_perimeter,'10000');
    target_value = local_B_at_point( ...
        kind{1},'0.999','0.04','10000');
    assert(I_inf(old_value) > I_sup(target_value));
    corrected_lower = compute_J_lower_bound( ...
        kind{1},'10000',area3148,perimeter3148);
    assert(corrected_lower <= I_inf(target_value));
end

% Cell 39970 has x_sup=1 as well, but its lower-left corner is rigorously
% above the floor.  It must retain the sharper lower-left comparison rule.
cell39970 = local_input_cell(input_table,39970);
assert(strcmp(cell39970.x_inf,'0.99902374999999999'));
assert(strcmp(cell39970.x_sup,'1'));
assert(strcmp(cell39970.theta_inf,'0.040905999999999998'));
assert(strcmp(cell39970.theta_sup,'0.042951249999999996'));
[area39970,perimeter39970,diagnostics39970] = ...
    local_geometry(cell39970);
assert(diagnostics39970.lower_corner_certified_above_floor);
assert(~diagnostics39970.use_floor_comparison);
assert(strcmp(diagnostics39970.comparison_geometry_rule, ...
    'lower-left parameter comparison point'));
local_assert_target_in_cell(cell39970,'0.9991','0.041');
local_assert_functional_comparison( ...
    area39970,perimeter39970,'0.9991','0.041');

% A cell lying wholly below the floor has no Omega_mid target.  The floor
% rule is still a safe (vacuous) comparison, and the two returned geometry
% evaluations need not be ordered by area in this case.
fully_below = struct( ...
    'x_inf','0.6','x_sup','0.61', ...
    'theta_inf','0.01','theta_sup','0.02');
[area_below,~,diagnostics_below] = local_geometry(fully_below);
assert(~diagnostics_below.lower_corner_certified_above_floor);
assert(diagnostics_below.use_floor_comparison);
assert(diagnostics_below.cell_certified_below_floor);
assert(I_inf(area_below(1)) > I_sup(area_below(2)));
local_assert_floor_comparison( ...
    area_below,diagnostics_below,fully_below.x_sup);

% Reproduce the production outward-overlap subdivision.  The lower-left
% descendants of cell 1363 become wholly below the floor at depth four;
% midpoint enclosures overlap, so no endpoint is lost between children.
node = cell1363;
for depth = 1:4
    children = local_split_like_production(node);
    local_assert_split_covers_parent(node,children);
    if depth == 1
        [~,~,upper_child] = local_geometry(children(3));
        assert(upper_child.lower_corner_certified_above_floor);
        assert(~upper_child.use_floor_comparison);
    end
    node = children(1);
end
[area_subdivision,perimeter_subdivision,subdivision_diagnostics] = ...
    local_geometry(node);
assert(subdivision_diagnostics.cell_certified_below_floor);
assert(subdivision_diagnostics.use_floor_comparison);
local_assert_floor_comparison( ...
    area_subdivision,subdivision_diagnostics,node.x_sup);
assert(all(isfinite([I_inf(perimeter_subdivision(1)), ...
    I_sup(perimeter_subdivision(1))])));

% A NaN must never select the lower-left branch through MATLAB's false
% comparison semantics.
nonfinite_rejected = false;
try
    compute_geometry_bounds(NaN,0.6,0.1,0.2);
catch ME
    nonfinite_rejected = strcmp( ...
        ME.identifier,'compute_geometry_bounds:NonfiniteGeometry');
end
assert(nonfinite_rejected);

fprintf('test_omega_mid_floor_geometry: PASS\n');
end


function [area,perimeter,diagnostics] = local_geometry(cell_data)
[area,perimeter,diagnostics] = compute_geometry_bounds( ...
    cell_data.x_inf,cell_data.x_sup, ...
    cell_data.theta_inf,cell_data.theta_sup);
end


function local_assert_floor_comparison(area,diagnostics,x_sup)
assert(~diagnostics.lower_corner_certified_above_floor);
assert(diagnostics.use_floor_comparison);
assert(diagnostics.bottom_straddle); % backward-compatible alias
assert(strcmp(diagnostics.comparison_geometry_rule, ...
    'right endpoint on y floor comparison point'));
assert(strcmp(diagnostics.lower_geometry_rule, ...
    diagnostics.comparison_geometry_rule));
x_sup_double = local_number(x_sup);
assert(I_inf(diagnostics.p_comparison.x) <= x_sup_double);
assert(I_sup(diagnostics.p_comparison.x) >= x_sup_double);
assert(I_inf(diagnostics.p_comparison.y) <= 0.04);
assert(I_sup(diagnostics.p_comparison.y) >= 0.04);
assert(I_inf(area(1)) <= 0.02 && I_sup(area(1)) >= 0.02);
assert(local_same_interval( ...
    diagnostics.p_lower.x,diagnostics.p_comparison.x));
assert(local_same_interval( ...
    diagnostics.p_lower.y,diagnostics.p_comparison.y));
end


function local_assert_target_in_cell( ...
    cell_data,x_value,y_value,strictly_above_floor)
if nargin < 4
    strictly_above_floor = true;
end
x = I_intval(x_value);
y = I_intval(y_value);
theta = atan(y/x);
x_lo = I_intval(cell_data.x_inf);
x_hi = I_intval(cell_data.x_sup);
t_lo = I_intval(cell_data.theta_inf);
t_hi = I_intval(cell_data.theta_sup);
assert(I_inf(x) > I_sup(x_lo));
assert(I_sup(x) < I_inf(x_hi));
assert(I_inf(theta) > I_sup(t_lo));
assert(I_sup(theta) < I_inf(t_hi));
if strictly_above_floor
    assert(I_inf(y) > I_sup(I_intval('0.04')));
else
    assert(local_same_interval(y,I_intval('0.04')));
end
assert(I_sup(x^2+y^2) < I_inf(I_intval('1')));
end


function local_assert_functional_comparison( ...
    area,perimeter,x_target,y_target)
for lambda = {'0','10000'}
    for kind = {'J1','J2'}
        [computed_lower,diagnostics] = compute_J_lower_bound( ...
            kind{1},lambda{1},area,perimeter);
        target_value = local_B_at_point( ...
            kind{1},x_target,y_target,lambda{1});
        assert(isfinite(computed_lower));
        assert(isfinite(I_inf(target_value)) ...
            && isfinite(I_sup(target_value)));
        assert(computed_lower <= I_inf(target_value));
        assert(strcmp(diagnostics.B_eval_point, ...
            'certified comparison point'));
    end
end
end


function value = local_B_at_point(kind,x_value,y_value,lambda)
x = I_intval(x_value);
y = I_intval(y_value);
area = y/I_intval('2');
one = I_intval('1');
perimeter = one+sqrt(x^2+y^2)+sqrt((one-x)^2+y^2);
value = local_B_from_geometry(kind,area,perimeter,lambda);
end


function value = local_B_from_geometry(kind,area,perimeter,lambda)
area = I_intval(area);
perimeter = I_intval(perimeter);
lambda = I_intval(lambda);
if strcmpi(kind,'J1')
    value = lambda*area ...
        -(I_pi^2/I_intval('16'))*perimeter^2/area ...
        -I_intval('7')*sqrt(I_intval('3'))*I_pi^2/I_intval('12');
elseif strcmpi(kind,'J2')
    constant = I_intval('4')*I_pi^2 ...
        /(I_intval('3')+sqrt(I_pi*sqrt(I_intval('3'))))^2;
    value = lambda*area ...
        -constant*(perimeter+sqrt(I_intval('4')*I_pi*area))^2 ...
        /(I_intval('4')*area);
else
    error('test_omega_mid_floor_geometry:UnknownFunctional', ...
        'Unknown functional %s.',kind);
end
end


function children = local_split_like_production(cell_data)
x_lo_interval = I_intval(cell_data.x_inf);
x_hi_interval = I_intval(cell_data.x_sup);
t_lo_interval = I_intval(cell_data.theta_inf);
t_hi_interval = I_intval(cell_data.theta_sup);
x_mid_interval = (x_lo_interval+x_hi_interval)/2;
t_mid_interval = (t_lo_interval+t_hi_interval)/2;
x_lo = I_inf(x_lo_interval);
x_hi = I_sup(x_hi_interval);
t_lo = I_inf(t_lo_interval);
t_hi = I_sup(t_hi_interval);
x_mid_lo = I_inf(x_mid_interval);
x_mid_hi = I_sup(x_mid_interval);
t_mid_lo = I_inf(t_mid_interval);
t_mid_hi = I_sup(t_mid_interval);

children = repmat(cell_data,4,1);
children(1).x_inf = x_lo;
children(1).x_sup = x_mid_hi;
children(1).theta_inf = t_lo;
children(1).theta_sup = t_mid_hi;
children(2).x_inf = x_mid_lo;
children(2).x_sup = x_hi;
children(2).theta_inf = t_lo;
children(2).theta_sup = t_mid_hi;
children(3).x_inf = x_lo;
children(3).x_sup = x_mid_hi;
children(3).theta_inf = t_mid_lo;
children(3).theta_sup = t_hi;
children(4).x_inf = x_mid_lo;
children(4).x_sup = x_hi;
children(4).theta_inf = t_mid_lo;
children(4).theta_sup = t_hi;
end


function local_assert_split_covers_parent(parent,children)
parent_x_lo = I_inf(I_intval(parent.x_inf));
parent_x_hi = I_sup(I_intval(parent.x_sup));
parent_t_lo = I_inf(I_intval(parent.theta_inf));
parent_t_hi = I_sup(I_intval(parent.theta_sup));
assert(children(1).x_inf == parent_x_lo);
assert(children(3).x_inf == parent_x_lo);
assert(children(2).x_sup == parent_x_hi);
assert(children(4).x_sup == parent_x_hi);
assert(children(1).theta_inf == parent_t_lo);
assert(children(2).theta_inf == parent_t_lo);
assert(children(3).theta_sup == parent_t_hi);
assert(children(4).theta_sup == parent_t_hi);
assert(children(1).x_sup >= children(2).x_inf);
assert(children(3).x_sup >= children(4).x_inf);
assert(children(1).theta_sup >= children(3).theta_inf);
assert(children(2).theta_sup >= children(4).theta_inf);
end


function table_value = local_read_input_table(input_file)
options = detectImportOptions(input_file,'TextType','string');
endpoints = {'x_inf','x_sup','theta_inf','theta_sup'};
options = setvartype(options,endpoints,'string');
table_value = readtable(input_file,options);
end


function cell_data = local_input_cell(table_value,cell_id)
raw_ids = table_value.i;
if isnumeric(raw_ids)
    ids = double(raw_ids);
else
    ids = str2double(string(raw_ids));
end
index = find(ids == cell_id);
assert(numel(index) == 1);
cell_data = struct( ...
    'x_inf',local_text(table_value.x_inf,index), ...
    'x_sup',local_text(table_value.x_sup,index), ...
    'theta_inf',local_text(table_value.theta_inf,index), ...
    'theta_sup',local_text(table_value.theta_sup,index));
end


function value = local_text(column,index)
raw = column(index);
if iscell(raw)
    raw = raw{1};
end
value = char(raw);
end


function value = local_number(input)
if ischar(input) || isstring(input)
    value = str2double(input);
else
    value = double(input);
end
end


function tf = local_same_interval(left,right)
tf = isequal(I_inf(left),I_inf(right)) ...
    && isequal(I_sup(left),I_sup(right));
end
