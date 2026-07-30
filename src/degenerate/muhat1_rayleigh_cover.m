function result = muhat1_rayleigh_cover(varargin)
%MUHAT1_RAYLEIGH_COVER  Verified Rayleigh cover for the thin-triangle model.
%
% For t>0 and s in [0,1], this routine bounds the first eigenvalue of
%
%   -u'' + V(t,s,x)u = muhat*u,  u in H^1_0(I^t),
%
% where I^t=[-t^(-2/3)(1+s),t^(-2/3)(1-s)] and V is the potential in
% Problem 2 of Endo--Liu (2025).  The change of variables
%
%   z=t^(2/3)x+s
%
% maps I^t onto [-1,1].  On each s-cell, a fixed test function
%
%   phi(z)=(1-z^2)*sum_j q_j P_{j-1}(z)
%
% is inserted into the Rayleigh quotient.  The vector q is computed at
% the cell midpoint and then frozen.  All integrals are exact polynomial
% antiderivatives; INTLAB is used only for outward-rounded evaluation.

opt = local_options(varargin{:});
is_interval = strcmp(opt.mode, 'interval');

if is_interval
    local_start_intlab();
    t = intval(opt.t_text);
    pi_value = intval('pi');
else
    t = str2double(opt.t_text);
    pi_value = pi;
end
t_double = str2double(opt.t_text);

basis = local_basis(opt.degree);
M_double = local_constant_matrix(basis.basis, basis.basis, -1, 1, false);
D_double = local_constant_matrix( ...
    basis.derivative, basis.derivative, -1, 1, false);
cells = repmat(local_empty_cell(), opt.s_cells, 1);
max_upper = -inf;
max_midpoint = -inf;
max_cell = 0;
max_down_scaled_upper = -inf;
max_down_cell = 0;

for cell_id = 1:opt.s_cells
    ds = (opt.s_hi-opt.s_lo)/opt.s_cells;
    s_lo = opt.s_lo+(cell_id-1)*ds;
    s_hi = opt.s_lo+cell_id*ds;
    s_mid = (s_lo+s_hi)/2;

    K_mid = local_stiffness( ...
        t_double, s_mid, basis, M_double, D_double, pi, false);
    [vectors, values] = eig(K_mid, M_double, 'vector');
    [mu_mid, order] = sort(real(values), 'ascend');
    q = real(vectors(:, order(1)));
    q = q/sqrt(q.'*M_double*q);
    if local_polyval(q.'*basis.legendre, 0) < 0
        q = -q;
    end

    if is_interval
        s_cell = infsup(s_lo, s_hi);
        rayleigh = local_scalar_rayleigh( ...
            t, s_cell, q, basis, pi_value, true);
        rayleigh_lower = inf(rayleigh);
        rayleigh_upper = sup(rayleigh);
        if ~isfinite(rayleigh_lower) || ~isfinite(rayleigh_upper)
            error('Nonfinite Rayleigh enclosure on s-cell %d: [%g,%g].', ...
                cell_id, s_lo, s_hi);
        end
    else
        rayleigh_lower = mu_mid(1);
        rayleigh_upper = mu_mid(1);
    end

    cells(cell_id).id = cell_id;
    cells(cell_id).s_lo = s_lo;
    cells(cell_id).s_hi = s_hi;
    cells(cell_id).s_mid = s_mid;
    cells(cell_id).midpoint_ritz = mu_mid(1);
    cells(cell_id).rayleigh_lower = rayleigh_lower;
    cells(cell_id).rayleigh_upper = rayleigh_upper;
    cells(cell_id).target_11_5_verified = is_interval && rayleigh_upper <= 11.5;
    if is_interval && ~isempty(opt.down_y_max)
        cells(cell_id).down_scaled_upper = local_down_scaled_upper( ...
            s_cell, rayleigh_upper, opt.down_y_max, ...
            opt.down_y_cells, pi_value);
        if cells(cell_id).down_scaled_upper > max_down_scaled_upper
            max_down_scaled_upper = cells(cell_id).down_scaled_upper;
            max_down_cell = cell_id;
        end
    end
    if opt.progress_every > 0 && mod(cell_id, opt.progress_every) == 0
        fprintf('completed %d/%d cells; current upper %.12g\n', ...
            cell_id, opt.s_cells, rayleigh_upper);
    end

    if rayleigh_upper > max_upper
        max_upper = rayleigh_upper;
        max_cell = cell_id;
    end
    max_midpoint = max(max_midpoint, mu_mid(1));
end

function rayleigh = local_scalar_rayleigh( ...
        t, s, q, basis, pi_value, is_interval)
if is_interval
    p = intval(q.')*intval(basis.legendre);
    phi = local_conv(p, intval([1 0 -1]));
    dphi = local_derivative(phi);
    p_squared = local_conv(p, p);
    left_integrand = local_conv(p_squared, intval([1 -2 1]));
    right_integrand = local_conv(p_squared, intval([1 2 1]));
    lo = intval(-1);
    hi = intval(1);
else
    p = q.'*basis.legendre;
    phi = conv(p, [1 0 -1]);
    dphi = local_derivative(phi);
    p_squared = conv(p, p);
    left_integrand = conv(p_squared, [1 -2 1]);
    right_integrand = conv(p_squared, [1 2 1]);
    lo = -1;
    hi = 1;
end

mass = local_integral(local_conv_or_conv(phi, phi, is_interval), lo, hi);
energy = local_integral( ...
    local_conv_or_conv(dphi, dphi, is_interval), lo, hi);
if is_interval
    left_primitive = local_primitive(left_integrand);
    right_primitive = local_primitive(right_integrand);
    left_polynomial = left_primitive;
    left_polynomial(1) = left_polynomial(1) ...
        -local_polyval(left_primitive, lo);
    right_polynomial = -right_primitive;
    right_polynomial(1) = right_polynomial(1) ...
        +local_polyval(right_primitive, hi);
    weighted_polynomial = local_poly_add( ...
        local_conv(left_polynomial, intval([1 2 1])), ...
        local_conv(right_polynomial, intval([1 -2 1])));
    weighted_polynomial(1) = weighted_polynomial(1)-mass;
    slope_polynomial = local_poly_add( ...
        left_polynomial, right_polynomial);
    a = t^(2/3);
    c = 3+4*pi_value^2;
    numerator_polynomial = ...
        (pi_value^2/a)*weighted_polynomial ...
        +(c*a^2/12)*local_pad( ...
            slope_polynomial, local_vector_length(weighted_polynomial));
    numerator_polynomial(1) = numerator_polynomial(1)+a^2*energy;
    numerator = local_polyval_centered(numerator_polynomial, s);
else
    left_value = local_integral(left_integrand, lo, s);
    right_value = local_integral(right_integrand, s, hi);
    a = t^(2/3);
    c = 3+4*pi_value^2;
    numerator = a^2*energy ...
        +(pi_value^2/a)*((1+s)^2*left_value ...
            +(1-s)^2*right_value-mass) ...
        +(c*a^2/12)*(left_value+right_value);
end
rayleigh = numerator/mass;
end

function product = local_conv_or_conv(left, right, is_interval)
if is_interval
    product = local_conv(left, right);
else
    product = conv(left, right);
end
end

result.mode = opt.mode;
result.rigor = local_if(is_interval, ...
    'certified_interval_rayleigh_upper', 'exploratory_double');
result.t = str2double(opt.t_text);
result.degree = opt.degree;
result.s_cells = opt.s_cells;
result.s_domain = [opt.s_lo, opt.s_hi];
result.max_rayleigh_upper = max_upper;
result.max_midpoint_ritz = max_midpoint;
result.max_cell_id = max_cell;
result.max_cell = [cells(max_cell).s_lo, cells(max_cell).s_hi];
result.target_11_5_verified = is_interval && max_upper <= 11.5;
if is_interval && ~isempty(opt.down_y_max)
    result.down_y_max = opt.down_y_max;
    result.down_y_cells = opt.down_y_cells;
    result.max_down_scaled_upper = max_down_scaled_upper;
    result.max_down_cell_id = max_down_cell;
    result.max_down_cell = [ ...
        cells(max_down_cell).s_lo, cells(max_down_cell).s_hi];
    result.down_conjecture_verified = max_down_scaled_upper < 0;
else
    result.down_y_max = [];
    result.down_y_cells = 0;
    result.max_down_scaled_upper = NaN;
    result.max_down_cell_id = 0;
    result.max_down_cell = [];
result.down_conjecture_verified = false;
end
result.cells = cells;

function upper = local_down_scaled_upper( ...
        s, rayleigh_upper, y_max, y_cells, pi_value)
one = intval(1);
mu = infsup(rayleigh_upper, rayleigh_upper);
y_limit = intval(sprintf('%.17g', y_max));
r_limit = sup(y_limit^(1/3));
nodes = linspace(0, r_limit, y_cells+1);
x = (one+s)/2;
r = infsup(nodes(1:end-1), nodes(2:end));
y = r.^3;
perimeter = one+sqrt(x^2+y.^2)+sqrt((one-x)^2+y.^2);
scaled = 6*(2*one)^(2/3)*mu*r.^2+6*pi_value^2 ...
    -2*pi_value^2*perimeter.^2 ...
    -4*sqrt(3*one)*pi_value^2*y;
upper = max(max(sup(scaled)));
end
end

function K = local_stiffness(t, s, basis, M, D, pi_value, is_interval)
a = t^(2/3);
n = numel(basis.legendre_cell);
if is_interval
    L = intval(zeros(n));
    R = intval(zeros(n));
else
    L = zeros(n);
    R = zeros(n);
end
for i = 1:n
    for j = i:n
        if is_interval
            product = local_conv( ...
                intval(basis.legendre_cell{i}), ...
                intval(basis.legendre_cell{j}));
            left_integrand = local_conv(product, intval([1 -2 1]));
            right_integrand = local_conv(product, intval([1 2 1]));
        else
            product = conv(basis.legendre_cell{i}, basis.legendre_cell{j});
            left_integrand = conv(product, [1 -2 1]);
            right_integrand = conv(product, [1 2 1]);
        end
        lij = local_integral(left_integrand, -1, s);
        rij = local_integral(right_integrand, s, 1);
        L(i,j) = lij;
        L(j,i) = lij;
        R(i,j) = rij;
        R(j,i) = rij;
    end
end
c = 3+4*pi_value^2;
K = a^2*D ...
    +(pi_value^2/a)*((1+s)^2*L+(1-s)^2*R-M) ...
    +(c*a^2/12)*(L+R);
end

function basis = local_basis(degree)
n = degree;
P = cell(n,1);
P{1} = 1;
if n >= 2
    P{2} = [0 1];
end
for k = 1:n-2
    zPk = [0 P{k+1}];
    Pkm1 = [P{k} zeros(1, numel(zPk)-numel(P{k}))];
    P{k+2} = ((2*k+1)*zPk-k*Pkm1)/(k+1);
end

b = cell(n,1);
db = cell(n,1);
max_length = n+2;
legendre = zeros(n,max_length);
for k = 1:n
    b{k} = conv(P{k}, [1 0 -1]);
    db{k} = local_derivative(b{k});
    legendre(k,1:numel(P{k})) = P{k};
end
basis.legendre_cell = P;
basis.legendre = legendre;
basis.basis = b;
basis.derivative = db;
end

function A = local_constant_matrix(left, right, lo, hi, is_interval)
n = numel(left);
if is_interval
    A = intval(zeros(n));
else
    A = zeros(n);
end
for i = 1:n
    for j = i:n
        if is_interval
            product = local_conv(intval(left{i}), intval(right{j}));
        else
            product = conv(left{i}, right{j});
        end
        value = local_integral(product, lo, hi);
        A(i,j) = value;
        A(j,i) = value;
    end
end
end

function product = local_conv(left, right)
n_left = local_vector_length(left);
n_right = local_vector_length(right);
product = left(1)*right(1)*zeros(1, n_left+n_right-1);
for i = 1:n_left
    for j = 1:n_right
        product(i+j-1) = product(i+j-1)+left(i)*right(j);
    end
end
end

function value = local_integral(coeff, lo, hi)
primitive = local_primitive(coeff);
value = local_polyval(primitive, hi)-local_polyval(primitive, lo);
end

function primitive = local_primitive(coeff)
n = local_vector_length(coeff);
primitive = [0 coeff./(1:n)];
end

function sum_value = local_poly_add(left, right)
n = max(local_vector_length(left), local_vector_length(right));
sum_value = local_pad(left, n)+local_pad(right, n);
end

function padded = local_pad(value, n)
n_value = local_vector_length(value);
padded = value(1)*zeros(1,n);
padded(1:n_value) = value;
end

function value = local_polyval(coeff, x)
n = local_vector_length(coeff);
value = coeff(n);
for k = n-1:-1:1
    value = value*x+coeff(k);
end
end

function value = local_polyval_centered(coeff, x)
n = local_vector_length(coeff);
x_mid = intval((inf(x)+sup(x))/2);
shifted = coeff(1)*zeros(1,n);
for j = 0:n-1
    for k = j:n-1
        shifted(j+1) = shifted(j+1) ...
            +coeff(k+1)*nchoosek(k,j)*x_mid^(k-j);
    end
end
value = local_polyval(shifted, x-x_mid);
end

function d = local_derivative(coeff)
n = local_vector_length(coeff);
if n == 1
    d = 0;
else
    d = (1:n-1).*coeff(2:n);
end
end

function n = local_vector_length(value)
dimensions = size(value);
n = max(dimensions);
end

function options = local_options(varargin)
p = inputParser;
addParameter(p, 'mode', 'double');
addParameter(p, 't', '0.38');
addParameter(p, 'degree', 14);
addParameter(p, 's_cells', 400);
addParameter(p, 'progress_every', 0);
addParameter(p, 's_lo', 0);
addParameter(p, 's_hi', 1);
addParameter(p, 'down_y_max', []);
addParameter(p, 'down_y_cells', 20);
parse(p, varargin{:});
options = p.Results;
options.mode = lower(char(options.mode));
if ~ismember(options.mode, {'double','interval'})
    error('mode must be ''double'' or ''interval''.');
end
if isnumeric(options.t)
    options.t_text = sprintf('%.17g', options.t);
else
    options.t_text = char(options.t);
end
if ~(str2double(options.t_text) > 0)
    error('t must be positive.');
end
if options.degree < 2 || options.degree ~= floor(options.degree)
    error('degree must be an integer at least 2.');
end
if options.s_cells < 1 || options.s_cells ~= floor(options.s_cells)
    error('s_cells must be a positive integer.');
end
if options.s_lo < 0 || options.s_hi > 1 || ...
        options.s_lo >= options.s_hi
    error('Require 0 <= s_lo < s_hi <= 1.');
end
if options.progress_every < 0 || ...
        options.progress_every ~= floor(options.progress_every)
    error('progress_every must be a nonnegative integer.');
end
if ~isempty(options.down_y_max) && options.down_y_max <= 0
    error('down_y_max must be positive.');
end
if options.down_y_cells < 1 || ...
        options.down_y_cells ~= floor(options.down_y_cells)
    error('down_y_cells must be a positive integer.');
end
end

function local_start_intlab()
if exist('intval', 'file') == 2
    return
end
intlab_root = getenv('INTLAB_ROOT');
if isempty(intlab_root)
    error('Set INTLAB_ROOT to the external INTLAB installation.');
end
addpath(intlab_root);
evalc('startintlab');
if exist('intval', 'file') ~= 2
    error('INTLAB initialization failed.');
end
end

function cell_result = local_empty_cell()
cell_result = struct( ...
    'id', 0, ...
    's_lo', NaN, ...
    's_hi', NaN, ...
    's_mid', NaN, ...
    'midpoint_ritz', NaN, ...
    'rayleigh_lower', NaN, ...
    'rayleigh_upper', NaN, ...
    'target_11_5_verified', false, ...
    'down_scaled_upper', NaN);
end

function value = local_if(condition, yes_value, no_value)
if condition
    value = yes_value;
else
    value = no_value;
end
end
