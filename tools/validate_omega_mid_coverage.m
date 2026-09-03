function certificate = validate_omega_mid_coverage(input_file,varargin)
%VALIDATE_OMEGA_MID_COVERAGE Rigorous geometric cover of Omega_mid.
%
% certificate = validate_omega_mid_coverage(input_file) proves that the
% closed (x,theta) rectangles in input_file cover
%
%   { (x,y): x >= 1/2, x^2+y^2 <= 1, y >= 0.04,
%              |(x,y)-(1/2,sqrt(3)/2)| >= 0.122 }.
%
% Covering the circle boundary is stronger than is needed if "outside the
% closed ball" is read as the strict complement.  CSV endpoint columns are
% imported as text.  Their exact decimal values determine the rectangular
% arrangement, and every analytic comparison is made with outward-rounded
% INTLAB enclosures constructed directly from those strings.  No geometric
% tolerance is used.
%
% Name/value option:
%   'manifest_file'  Write the returned certificate atomically as JSON.

% The algorithm uses the exact decimal coordinate arrangement.  A 2-D
% difference array determines which elementary (x,theta) rectangles are in
% the input union.  On each elementary x slab the active theta components
% are constant.  The required lower boundary atan(0.04/x) is decreasing.
% The upper boundary is the minimum of acos(x) and the lower branch of the
% excluded circle.  Its only possible interior extrema occur at the unique
% circle-branch critical point and at the circle/unit-circle switch; these
% algebraic points are included with interval enclosures below.

p = inputParser;
p.CaseSensitive = true;
if isprop(p,'PartialMatching')
    p.PartialMatching = false;
end
addRequired(p,'input_file',@(x) ischar(x) || isstring(x));
addParameter(p,'manifest_file','', ...
    @(x) ischar(x) || isstring(x));
parse(p,input_file,varargin{:});
input_file = char(p.Results.input_file);
manifest_file = char(p.Results.manifest_file);

global INTERVAL_MODE
if isempty(INTERVAL_MODE) || ~logical(INTERVAL_MODE)
    error('validate_omega_mid_coverage:IntervalModeRequired', ...
        ['A geometric proof certificate requires initialized INTLAB ', ...
         'interval mode.']);
end
probe = I_intval('0.5');
if isnumeric(probe)
    error('validate_omega_mid_coverage:IntervalConstructionFailed', ...
        'I_intval returned a non-interval value in interval mode.');
end
if ~exist(input_file,'file')
    error('validate_omega_mid_coverage:MissingInput', ...
        'The cell-definition CSV does not exist: %s',input_file);
end

input_sha256 = local_file_sha256(input_file);
table_value = local_read_endpoint_table(input_file);
cell_count = height(table_value);
if cell_count == 0
    error('validate_omega_mid_coverage:EmptyInput', ...
        'The cell-definition CSV is empty.');
end

[x_text,x_lo_index,x_hi_index] = local_coordinate_axis( ...
    table_value.x_inf,table_value.x_sup,'x');
[theta_text,theta_lo_index,theta_hi_index] = ...
    local_coordinate_axis( ...
        table_value.theta_inf,table_value.theta_sup,'theta');
if any(x_lo_index >= x_hi_index) ...
        || any(theta_lo_index >= theta_hi_index)
    error('validate_omega_mid_coverage:NonpositiveCellWidth', ...
        ['Every cell must have strictly positive x and theta width in ', ...
         'the exact decimal ordering.']);
end

[x_interval,x_lower,x_upper] = ...
    local_interval_coordinates(x_text,'x');
[~,theta_lower,theta_upper] = ...
    local_interval_coordinates(theta_text,'theta');

zero = I_intval('0');
half = I_intval('0.5');
one = I_intval('1');
three = I_intval('3');
y_floor = I_intval('0.04');
radius = I_intval('0.122');
p0_y = sqrt(three)/I_intval('2');
theta_cap = acos(half);
x_floor = sqrt(one-y_floor^2);
x_ball_right = half+radius;

% Exact cells may extend to the curved Omega boundary in theta, but their
% rectangular parameter endpoints must remain in the fundamental box.
if x_lower(1) < local_lower(half) ...
        || x_upper(end) > local_upper(one) ...
        || theta_lower(1) < local_lower(zero) ...
        || theta_upper(end) > local_upper(theta_cap)
    error('validate_omega_mid_coverage:CellOutsideParameterBox', ...
        ['Cell endpoints must lie in [1/2,1] x [0,pi/3]. ', ...
         'Extra rectangles outside that box are not accepted.']);
end
if x_upper(1) > local_lower(half) ...
        || x_lower(end) < local_upper(x_floor)
    error('validate_omega_mid_coverage:XDomainNotEnclosed', ...
        ['The cell coordinate range does not enclose the complete target ', ...
         'x interval [1/2,sqrt(1-0.04^2)].']);
end

% Algebraic break points of the upper theta boundary.
a = I_intval('2')-radius^2;
x_switch = (a+sqrt(three*(I_intval('4')-a^2)))/I_intval('4');
d_critical = I_intval('0.5')*( ...
    -radius^2+radius*sqrt(three*(one-radius^2)));
x_ball_critical = half+d_critical;
if ~(local_upper(y_floor) < local_lower(p0_y-radius) ...
        && local_upper(half) < local_lower(x_ball_critical) ...
        && local_upper(x_ball_critical) < local_lower(x_switch) ...
        && local_upper(x_switch) < local_lower(x_ball_right) ...
        && local_upper(x_ball_right) < local_lower(x_floor))
    error('validate_omega_mid_coverage:UnresolvedBoundaryOrdering', ...
        ['INTLAB did not certify the required ordering of the circle ', ...
         'critical point, boundary switch, and target endpoint.']);
end

nx = numel(x_text);
nt = numel(theta_text);
event_row = [x_lo_index;x_hi_index;x_lo_index;x_hi_index];
event_column = [theta_lo_index;theta_lo_index; ...
    theta_hi_index;theta_hi_index];
event_value = [ones(cell_count,1);-ones(cell_count,1); ...
    -ones(cell_count,1);ones(cell_count,1)];
difference = accumarray( ...
    [event_row,event_column],event_value,[nx,nt],@sum,0,false);
cover_count = cumsum(cumsum(difference,1),2);
if any(cover_count(:) < 0) ...
        || any(cover_count(:) ~= floor(cover_count(:)))
    error('validate_omega_mid_coverage:InvalidArrangementCount', ...
        'The rectangle-arrangement coverage counts are inconsistent.');
end
covered = cover_count(1:nx-1,1:nt-1) > 0;
clear difference cover_count event_row event_column event_value

critical_points = {x_ball_critical,x_switch,x_ball_right};
complete = true;
witness = struct([]);
checked_x_slabs = 0;
maximum_components = 0;
minimum_lower_margin = Inf;
minimum_upper_margin = Inf;
candidate_evaluations = 0;

for ix = 1:nx-1
    x_left = x_interval{ix};
    x_right_csv = x_interval{ix+1};

    % There is no target interior to the left of x=1/2 or to the right of
    % the point where the unit-circle height equals 0.04.
    if local_upper(x_right_csv) <= local_lower(half) ...
            || local_lower(x_left) >= local_upper(x_floor)
        continue;
    end
    checked_x_slabs = checked_x_slabs+1;
    if local_upper(x_right_csv) <= local_lower(x_floor)
        x_right = x_right_csv;
    else
        % This also safely handles an interval-indeterminate comparison at
        % the irrational endpoint: evaluating as far as x_floor can only
        % enlarge the required enclosure by roundoff width.
        x_right = x_floor;
    end

    lower_theta = atan(y_floor/x_right);
    required_lower_floor = local_lower(lower_theta);

    required_upper_ceiling = max( ...
        local_upper_boundary_at(x_left,half,radius,p0_y, ...
            x_ball_right), ...
        local_upper_boundary_at(x_right,half,radius,p0_y, ...
            x_ball_right));
    candidate_evaluations = candidate_evaluations+2;
    for q = 1:numel(critical_points)
        candidate = critical_points{q};
        if local_intervals_overlap( ...
                candidate,x_left,x_right)
            required_upper_ceiling = max(required_upper_ceiling, ...
                local_upper_boundary_at( ...
                    candidate,half,radius,p0_y,x_ball_right));
            candidate_evaluations = candidate_evaluations+1;
        end
    end

    row_covered = covered(ix,:);
    component_start = find(row_covered ...
        & [true,~row_covered(1:end-1)]);
    component_end = find(row_covered ...
        & [~row_covered(2:end),true])+1;
    maximum_components = max(maximum_components,numel(component_start));

    slab_ok = false;
    best_lower_margin = -Inf;
    best_upper_margin = -Inf;
    for q = 1:numel(component_start)
        lower_margin = required_lower_floor ...
            -theta_upper(component_start(q));
        upper_margin = theta_lower(component_end(q)) ...
            -required_upper_ceiling;
        if lower_margin >= 0 && upper_margin >= 0
            slab_ok = true;
            minimum_lower_margin = min( ...
                minimum_lower_margin,lower_margin);
            minimum_upper_margin = min( ...
                minimum_upper_margin,upper_margin);
            break;
        end
        best_lower_margin = max(best_lower_margin,lower_margin);
        best_upper_margin = max(best_upper_margin,upper_margin);
    end

    if ~slab_ok
        complete = false;
        witness = struct( ...
            'kind','uncovered_elementary_x_slab', ...
            'x_inf_decimal',x_text{ix}, ...
            'x_sup_decimal',x_text{ix+1}, ...
            'x_slab_ordinal',ix, ...
            'required_theta_lower_certified_floor', ...
                required_lower_floor, ...
            'required_theta_upper_certified_ceiling', ...
                required_upper_ceiling, ...
            'covered_theta_component_count', ...
                numel(component_start), ...
            'best_lower_margin',best_lower_margin, ...
            'best_upper_margin',best_upper_margin);
        break;
    end
end

if checked_x_slabs == 0
    error('validate_omega_mid_coverage:NoTargetSlabsChecked', ...
        'No elementary x slab intersected the target region.');
end
if ~complete
    minimum_lower_margin = NaN;
    minimum_upper_margin = NaN;
end
if ~strcmp(local_file_sha256(input_file),input_sha256)
    error('validate_omega_mid_coverage:InputChangedDuringProof', ...
        ['The CSV changed after its initial hash and before certificate ', ...
         'publication; no mixed-state coverage proof was emitted.']);
end

certificate = struct();
certificate.schema = ...
    'lowerboundsineq.omega_mid_geometric_coverage.v1';
certificate.algorithm = ...
    'exact-decimal-arrangement-intlab-boundary-extrema-v1';
certificate.status = local_ternary(complete,'certified','failed');
certificate.complete = logical(complete);
certificate.rigorous = true;
certificate.input = struct( ...
    'sha256',input_sha256, ...
    'cell_count',cell_count, ...
    'endpoint_import','CSV text', ...
    'interval_conversion','decimal text directly to INTLAB');
certificate.region = struct( ...
    'omega','x>=1/2, y>0, x^2+y^2<=1', ...
    'y_floor_decimal','0.04', ...
    'ball_center_x_decimal','0.5', ...
    'ball_center_y','sqrt(3)/2', ...
    'ball_radius_decimal','0.122', ...
    'certified_distance_relation','distance >= 0.122', ...
    'requested_strict_complement_also_covered',true, ...
    'coordinate_map','theta=atan(y/x)');
certificate.proof = struct( ...
    'decimal_topology_exact',true, ...
    'binary_double_endpoint_roundtrip',false, ...
    'arbitrary_tolerance_used',false, ...
    'closed_rectangles',true, ...
    'elementary_grid_union_checked',true, ...
    'lower_boundary_monotonicity','atan(0.04/x) decreasing', ...
    'upper_boundary_extrema', ...
        ['slab endpoints plus the enclosed circle-branch critical ', ...
         'point, circle/unit-circle switch, and circle right endpoint'], ...
    'interval_comparison_rule', ...
        ['cell lower endpoint upper bound <= required lower floor; ', ...
         'cell upper endpoint lower bound >= required upper ceiling']);
certificate.statistics = struct( ...
    'unique_x_coordinates',nx, ...
    'unique_theta_coordinates',nt, ...
    'elementary_x_slabs',nx-1, ...
    'target_x_slabs_checked',checked_x_slabs, ...
    'maximum_covered_theta_components_per_slab',maximum_components, ...
    'upper_boundary_candidate_evaluations',candidate_evaluations, ...
    'minimum_certified_lower_margin',minimum_lower_margin, ...
    'minimum_certified_upper_margin',minimum_upper_margin);
certificate.witness = witness;

witness_key = 'none';
if ~complete
    witness_key = sprintf('%s,%s,%d', ...
        witness.x_inf_decimal,witness.x_sup_decimal, ...
        witness.x_slab_ordinal);
end
payload = sprintf([ ...
    'schema=%s|algorithm=%s|input=%s|cells=%d|', ...
    'region=Omega,x>=0.5,x^2+y^2<=1,y>=0.04,', ...
    'distance>=0.122|xcoords=%d|tcoords=%d|slabs=%d|', ...
    'complete=%d|witness=%s'], ...
    certificate.schema,certificate.algorithm,input_sha256,cell_count, ...
    nx,nt,checked_x_slabs,complete,witness_key);
certificate.hashes = struct( ...
    'algorithm','SHA-256', ...
    'input_csv',input_sha256, ...
    'proof_payload',local_sha256_utf8(payload));

if ~isempty(manifest_file)
    local_write_json_atomic(manifest_file,certificate);
end
end


function table_value = local_read_endpoint_table(input_file)
import_options = detectImportOptions(input_file,'TextType','string');
required = {'x_inf','x_sup','theta_inf','theta_sup'};
missing = required(~ismember(required,import_options.VariableNames));
if ~isempty(missing)
    error('validate_omega_mid_coverage:MissingEndpointColumn', ...
        'Missing endpoint column %s in %s.',missing{1},input_file);
end
import_options = setvartype(import_options,required,'string');
table_value = readtable(input_file,import_options);
end


function [coordinates,lo_index,hi_index] = local_coordinate_axis( ...
    lower_column,upper_column,axis_name)
lower_text = local_column_text(lower_column,axis_name);
upper_text = local_column_text(upper_column,axis_name);
raw = [lower_text;upper_text];
[raw_unique,~,raw_map] = unique(raw);
canonical_raw = cell(size(raw_unique));
for k = 1:numel(raw_unique)
    canonical_raw{k} = local_canonical_decimal( ...
        raw_unique{k},axis_name);
end
[canonical_unique,~,canonical_map] = unique(canonical_raw);
value_map = canonical_map(raw_map);
order = local_decimal_sort_order(canonical_unique);
coordinates = canonical_unique(order);
inverse_order = zeros(numel(order),1);
inverse_order(order) = 1:numel(order);
sorted_map = inverse_order(value_map);
n = numel(lower_text);
lo_index = sorted_map(1:n);
hi_index = sorted_map(n+1:end);
end


function values = local_column_text(column,axis_name)
if isstring(column)
    values = cellstr(column(:));
elseif iscell(column)
    values = cell(numel(column),1);
    for k = 1:numel(column)
        if ~(ischar(column{k}) || isstring(column{k}))
            error('validate_omega_mid_coverage:EndpointNotText', ...
                '%s endpoints were not imported as text.',axis_name);
        end
        values{k} = char(column{k});
    end
elseif ischar(column)
    values = cellstr(column);
else
    error('validate_omega_mid_coverage:EndpointNotText', ...
        '%s endpoints were not imported as text.',axis_name);
end
values = values(:);
end


function canonical = local_canonical_decimal(value,axis_name)
value = strtrim(char(value));
if isempty(value) || isempty(regexp(value, ...
        '^[+]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)$','once'))
    error('validate_omega_mid_coverage:InvalidDecimalEndpoint', ...
        ['%s endpoint "%s" is not a finite plain decimal.  Scientific ', ...
         'notation is rejected so the exact decimal topology is explicit.'], ...
        axis_name,value);
end
if value(1) == '+'
    value = value(2:end);
end
dot = find(value == '.',1);
if isempty(dot)
    integer = value;
    fraction = '';
else
    integer = value(1:dot-1);
    fraction = value(dot+1:end);
end
if isempty(integer)
    integer = '0';
end
first_nonzero = find(integer ~= '0',1);
if isempty(first_nonzero)
    integer = '0';
else
    integer = integer(first_nonzero:end);
end
while ~isempty(fraction) && fraction(end) == '0'
    fraction(end) = [];
end
if isempty(fraction)
    canonical = integer;
else
    canonical = [integer,'.',fraction];
end
end


function order = local_decimal_sort_order(values)
n = numel(values);
integer = cell(n,1);
fraction = cell(n,1);
maximum_integer = 1;
maximum_fraction = 0;
for k = 1:n
    dot = find(values{k} == '.',1);
    if isempty(dot)
        integer{k} = values{k};
        fraction{k} = '';
    else
        integer{k} = values{k}(1:dot-1);
        fraction{k} = values{k}(dot+1:end);
    end
    maximum_integer = max(maximum_integer,length(integer{k}));
    maximum_fraction = max(maximum_fraction,length(fraction{k}));
end
keys = cell(n,1);
for k = 1:n
    keys{k} = [ ...
        repmat('0',1,maximum_integer-length(integer{k})),integer{k},'.', ...
        fraction{k}, ...
        repmat('0',1,maximum_fraction-length(fraction{k}))];
end
[~,order] = sortrows(char(keys));
end


function [intervals,lower,upper] = ...
    local_interval_coordinates(decimal_text,axis_name)
n = numel(decimal_text);
intervals = cell(n,1);
lower = zeros(n,1);
upper = zeros(n,1);
for k = 1:n
    % This is the only decimal-to-number conversion for endpoint geometry.
    % In particular, there is deliberately no str2double or numeric table
    % conversion before INTLAB sees the decimal string.
    intervals{k} = I_intval(decimal_text{k});
    if isnumeric(intervals{k})
        error('validate_omega_mid_coverage:IntervalConstructionFailed', ...
            '%s endpoint did not produce an INTLAB interval.',axis_name);
    end
    lower(k) = local_lower(intervals{k});
    upper(k) = local_upper(intervals{k});
    if ~isfinite(lower(k)) || ~isfinite(upper(k)) ...
            || lower(k) > upper(k)
        error('validate_omega_mid_coverage:NonfiniteIntervalEndpoint', ...
            '%s endpoint did not have a finite interval enclosure.', ...
            axis_name);
    end
end
end


function bound = local_upper_boundary_at( ...
    x,half,radius,p0_y,x_ball_right)
omega_theta = acos(x);
bound = local_upper(omega_theta);

% For x at or to the right of 1/2+r the unit-circle cap alone is used.
% A strict enclosure test avoids taking sqrt of a radicand whose lower
% interval endpoint can be negative solely at the exact circle endpoint.
if local_upper(x) < local_lower(x_ball_right)
    radicand = radius^2-(x-half)^2;
    if local_lower(radicand) < 0
        error('validate_omega_mid_coverage:UnresolvedCircleRadicand', ...
            'The excluded-circle radicand was not certified nonnegative.');
    end
    ball_theta = atan((p0_y-sqrt(radicand))/x);
    bound = min(bound,local_upper(ball_theta));
end
if ~isfinite(bound)
    error('validate_omega_mid_coverage:NonfiniteBoundaryBound', ...
        'An upper-boundary interval was nonfinite.');
end
end


function tf = local_intervals_overlap(candidate,left,right)
tf = ~(local_upper(candidate) < local_lower(left) ...
    || local_lower(candidate) > local_upper(right));
end


function value = local_lower(interval_value)
value = double(inf(interval_value));
end


function value = local_upper(interval_value)
value = double(sup(interval_value));
end


function digest = local_file_sha256(filename)
fid = fopen(filename,'rb');
if fid < 0
    error('validate_omega_mid_coverage:HashOpenFailed', ...
        'Cannot open %s for hashing.',filename);
end
cleanup = onCleanup(@() fclose(fid));
engine = javaMethod( ...
    'getInstance','java.security.MessageDigest','SHA-256');
while true
    bytes = fread(fid,1024*1024,'*uint8');
    if isempty(bytes)
        break;
    end
    engine.update(typecast(bytes(:),'int8'));
end
raw = typecast(engine.digest(),'uint8');
digest = lower(reshape(dec2hex(raw,2).',1,[]));
clear cleanup
end


function digest = local_sha256_utf8(value)
try
    bytes = unicode2native(char(value),'UTF-8');
catch
    bytes = uint8(char(value));
end
engine = javaMethod( ...
    'getInstance','java.security.MessageDigest','SHA-256');
bytes = uint8(bytes(:));
if ~isempty(bytes)
    engine.update(typecast(bytes,'int8'));
end
raw = typecast(engine.digest(),'uint8');
digest = lower(reshape(dec2hex(raw,2).',1,[]));
end


function local_write_json_atomic(filename,value)
directory = fileparts(filename);
if ~isempty(directory) && ~exist(directory,'dir')
    mkdir(directory);
end
try
    encoded = jsonencode(value,'PrettyPrint',true);
catch
    encoded = jsonencode(value);
end
if isempty(directory)
    directory = pwd;
end
temporary = [tempname(directory),'.json'];
cleanup = onCleanup(@() local_delete_if_present(temporary));
fid = fopen(temporary,'w');
if fid < 0
    error('validate_omega_mid_coverage:ManifestWriteFailed', ...
        'Cannot open a temporary coverage manifest.');
end
file_cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',encoded);
clear file_cleanup
[ok,message] = movefile(temporary,filename,'f');
if ~ok
    error('validate_omega_mid_coverage:ManifestPublishFailed', ...
        'Could not publish %s: %s',filename,message);
end
clear cleanup
end


function local_delete_if_present(filename)
if exist(filename,'file')
    delete(filename);
end
end


function value = local_ternary(condition,left,right)
if condition
    value = left;
else
    value = right;
end
end
