function data = triangle_bernstein_bubble_data(degree)
%TRIANGLE_BERNSTEIN_BUBBLE_DATA  Exact Bernstein data on the unit triangle.
%
% The trial space is
%
%   V_p = span{B^p_alpha : alpha_1,alpha_2,alpha_3 >= 1},
%
% where |alpha|=p and B^p_alpha is the barycentric Bernstein polynomial on
% conv{(0,0),(1,0),(0,1)}.  Hence V_p is a subspace of
% H^2(T) cap H^1_0(T).  The returned matrices contain only exact integer
% derivative coefficients and rigorously rounded rational Bernstein
% integrals.  They are independent of the physical triangle.

global INTERVAL_MODE
if isempty(INTERVAL_MODE)
    INTERVAL_MODE = 0;
end
if ~isscalar(degree) || ~isfinite(degree) ...
        || degree ~= floor(degree) || degree < 3
    error('triangle_bernstein_bubble_data:BadDegree', ...
        'degree must be an integer at least three.');
end

persistent cached_degree cached_interval_mode cached_data
mode = logical(INTERVAL_MODE);
if ~isempty(cached_data) && isequal(cached_degree,degree) ...
        && isequal(cached_interval_mode,mode)
    data = cached_data;
    return;
end

p = degree;
alpha_p = local_multiindices(p,false);
alpha_i = local_multiindices(p,true);
alpha_p1 = local_multiindices(p-1,false);
alpha_p2 = local_multiindices(p-2,false);

map_p = local_index_map(alpha_p,p);
map_p1 = local_index_map(alpha_p1,p-1);
map_p2 = local_index_map(alpha_p2,p-2);

n = size(alpha_i,1);
np = size(alpha_p,1);
np1 = size(alpha_p1,1);
np2 = size(alpha_p2,1);

E = I_zeros(np,n);
for col = 1:n
    E(map_p(local_key(alpha_i(col,:),p)),col) = I_intval(1);
end

Dx = I_zeros(np1,n);
Dy = I_zeros(np1,n);
for col = 1:n
    a = alpha_i(col,:);
    Dx = local_add_derivative_column( ...
        Dx,col,a,p,2,1,map_p1,p-1);
    Dy = local_add_derivative_column( ...
        Dy,col,a,p,3,1,map_p1,p-1);
end

Dxx = local_derivative_matrix_full(Dx,alpha_p1,p-1,2,1, ...
    map_p2,p-2);
Dxy = local_derivative_matrix_full(Dx,alpha_p1,p-1,3,1, ...
    map_p2,p-2);
Dyy = local_derivative_matrix_full(Dy,alpha_p1,p-1,3,1, ...
    map_p2,p-2);

Gpp = local_bernstein_gram(alpha_p,p,alpha_p,p);
Gp1p1 = local_bernstein_gram(alpha_p1,p-1,alpha_p1,p-1);
Gp2p2 = local_bernstein_gram(alpha_p2,p-2,alpha_p2,p-2);
Gp2p = local_bernstein_gram(alpha_p2,p-2,alpha_p,p);

M = E'*Gpp*E;
Axx = Dx'*Gp1p1*Dx;
Axy = Dx'*Gp1p1*Dy;
Ayy = Dy'*Gp1p1*Dy;

data = struct();
data.degree = p;
data.dimension = n;
data.alpha_interior = alpha_i;
data.alpha_full = alpha_p;
data.alpha_hessian = alpha_p2;
data.E = E;
data.Dx = Dx;
data.Dy = Dy;
data.Dxx = Dxx;
data.Dxy = Dxy;
data.Dyy = Dyy;
data.Gpp = Gpp;
data.Gp2p2 = Gp2p2;
data.Gp2p = Gp2p;
data.M = M;
data.Axx = Axx;
data.Axy = Axy;
data.Ayy = Ayy;

cached_degree = degree;
cached_interval_mode = mode;
cached_data = data;
end


function alpha = local_multiindices(p,interior_only)
rows = zeros((p+1)*(p+2)/2,3);
count = 0;
for a1 = 0:p
    for a2 = 0:(p-a1)
        a3 = p-a1-a2;
        if interior_only && any([a1,a2,a3] < 1)
            continue;
        end
        count = count+1;
        rows(count,:) = [a1,a2,a3];
    end
end
alpha = rows(1:count,:);
end


function map = local_index_map(alpha,p)
map = zeros((p+1)^3,1);
for row = 1:size(alpha,1)
    map(local_key(alpha(row,:),p)) = row;
end
end


function key = local_key(alpha,p)
key = 1+alpha(1)+(p+1)*alpha(2)+(p+1)^2*alpha(3);
end


function D = local_add_derivative_column( ...
    D,col,alpha,p,positive_index,negative_index,map_target,p_target)
% d/dx = p(B_{alpha-e_2}^{p-1}-B_{alpha-e_1}^{p-1}),
% d/dy = p(B_{alpha-e_3}^{p-1}-B_{alpha-e_1}^{p-1}).
if alpha(positive_index) > 0
    beta = alpha;
    beta(positive_index) = beta(positive_index)-1;
    row = map_target(local_key(beta,p_target));
    D(row,col) = D(row,col)+I_intval(p);
end
if alpha(negative_index) > 0
    beta = alpha;
    beta(negative_index) = beta(negative_index)-1;
    row = map_target(local_key(beta,p_target));
    D(row,col) = D(row,col)-I_intval(p);
end
end


function Dout = local_derivative_matrix_full( ...
    Din,alpha_source,p,positive_index,negative_index,map_target,p_target)
nsource = size(alpha_source,1);
Dfull = I_zeros(max(map_target),nsource);
for col = 1:nsource
    Dfull = local_add_derivative_column( ...
        Dfull,col,alpha_source(col,:),p, ...
        positive_index,negative_index,map_target,p_target);
end
Dout = Dfull*Din;
end


function G = local_bernstein_gram(alpha_p,p,alpha_q,q)
np = size(alpha_p,1);
nq = size(alpha_q,1);
G = I_zeros(np,nq);
denom_integral = (p+q+1)*(p+q+2);
for row = 1:np
    cp = local_multinomial(p,alpha_p(row,:));
    for col = 1:nq
        cq = local_multinomial(q,alpha_q(col,:));
        gamma = alpha_p(row,:)+alpha_q(col,:);
        cpq = local_multinomial(p+q,gamma);
        % Avoid forming cp*cq in floating point.  Each integer is exactly
        % representable at the degrees used here; interval multiplication
        % and division then enclose the exact rational integral.
        G(row,col) = I_intval(cp)*I_intval(cq) ...
            /I_intval(cpq)/I_intval(denom_integral);
    end
end
end


function value = local_multinomial(n,alpha)
value = 1;
remaining = n;
for k = 1:(numel(alpha)-1)
    value = value*nchoosek(remaining,alpha(k));
    remaining = remaining-alpha(k);
end
if value > flintmax
    error('triangle_bernstein_bubble_data:IntegerOverflow', ...
        'Increase the exact-integer implementation for this degree.');
end
end
