function certificate = triangle_spectral_certificate_cr_lg( ...
    triangle,N_LG,N_rho,bound_order,index)
%TRIANGLE_SPECTRAL_CERTIFICATE_CR_LG  Three verified spectral endpoints.
%
% For a point triangle T=conv{(0,0),(1,0),(x,y)}, return, for n=INDEX,
%
%   Lk_LG <= lambda_k(T) <= Uk_Ritz, 1<=k<=n,
%   Un_Ritz < L(n+1)_CR <= lambda_(n+1)(T).
%
% The final CR--Liu endpoint is the Lehmann--Goerisch shift.  No upper
% cluster dimension and no enclosure of lambda_(n+2) is required.

global INTERVAL_MODE
if isempty(INTERVAL_MODE)
    INTERVAL_MODE = 0;
end
if nargin < 5 || isempty(index)
    index = 1;
end
if ~isscalar(index) || ~isfinite(index) ...
        || index < 1 || index ~= floor(index)
    error('triangle_spectral_certificate_cr_lg:BadIndex', ...
        'INDEX must be a positive integer.');
end

triangle = I_intval(triangle);
if length(triangle) ~= 6
    error('triangle_spectral_certificate_cr_lg:BadTriangle', ...
        'TRIANGLE must have length six.');
end
prefix = [0,0,1,0];
for k = 1:4
    if I_inf(triangle(k)) ~= prefix(k) ...
            || I_sup(triangle(k)) ~= prefix(k)
        error('triangle_spectral_certificate_cr_lg:BadTriangle', ...
            'TRIANGLE must be conv{(0,0),(1,0),(x,y)}.');
    end
end
if I_inf(triangle(5)) ~= I_sup(triangle(5)) ...
        || I_inf(triangle(6)) ~= I_sup(triangle(6)) ...
        || ~(I_inf(triangle(6)) > 0)
    error('triangle_spectral_certificate_cr_lg:BadTriangle', ...
        'The reference geometry must be a point triangle of positive height.');
end

parameters = [N_LG,N_rho,bound_order];
if any(~isfinite(parameters)) || any(parameters < 1) ...
        || any(parameters ~= floor(parameters))
    error('triangle_spectral_certificate_cr_lg:BadParameters', ...
        'N_LG, N_rho, and BOUND_ORDER must be positive integers.');
end

[lambda_bounds,~,~,~,~,~,~,~,~,details] = ...
    calc_eigen_bounds_any_order_1k_wh( ...
        index,triangle,N_LG,N_rho,bound_order);
LG_lower = I_intval(zeros(index,1));
Ritz_upper = I_intval(zeros(index,1));
for k = 1:index
    LG_lower(k) = I_intval(I_inf(lambda_bounds(k)));
    Ritz_upper(k) = I_intval(I_sup(lambda_bounds(k)));
end
Lnext = I_intval(I_inf(details.CR_lambda_next_lower));
if any(I_inf(LG_lower) <= 0) ...
        || any(I_inf(LG_lower) > I_sup(Ritz_upper)) ...
        || ~(I_inf(Lnext) > I_sup(Ritz_upper(index)))
    error('triangle_spectral_certificate_cr_lg:InvalidEndpoints', ...
        ['Could not certify positive LG/Ritz enclosures and ', ...
         'U_i^{Ritz}<L_{i+1}^{CR}.']);
end

certificate = struct();
certificate.schema = ...
    'lowerboundsineq.triangle-spectral-cr-lg-certificate.v1';
certificate.method = ...
    'Lehmann-Goerisch+conforming-Ritz+CR-Liu';
certificate.reference_triangle = triangle;
certificate.index = index;
certificate.lambda_LG_lower = LG_lower;
certificate.lambda_Ritz_upper = Ritz_upper;
certificate.lambda_next_CR_lower = Lnext;
if index == 1
    certificate.lambda1_LG_lower = LG_lower(1);
    certificate.lambda1_Ritz_upper = Ritz_upper(1);
    certificate.lambda2_CR_lower = Lnext;
end
certificate.rigorous = logical(INTERVAL_MODE);
certificate.mesh_used = true;
certificate.N_LG = N_LG;
certificate.N_rho = N_rho;
certificate.bound_order = bound_order;
certificate.details = details;
end
