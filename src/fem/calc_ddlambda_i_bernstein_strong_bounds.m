function [lambda,dlambda,ddlambda_lower,ddlambda_upper,diagnostics] = ...
    calc_ddlambda_i_bernstein_strong_bounds( ...
        index,triangle,e_direction,degree,spectral_options)
%CALC_DDLAMBDA_I_BERNSTEIN_STRONG_BOUNDS  Verified lambda_i Hessian bound.
%
% INDEX must be a positive integer whose two adjacent spectral gaps can be
% certified.  The implementation uses global Bernstein strong residuals;
% only the scalar LG/Ritz/CR--Liu endpoints use a finite element mesh.

if nargin < 4 || isempty(degree)
    degree = 11;
end
if nargin < 5 || isempty(spectral_options)
    spectral_options = struct();
end
if ~isscalar(index) || ~isfinite(index) ...
        || index < 1 || index ~= floor(index)
    error('calc_ddlambda_i_bernstein_strong_bounds:BadIndex', ...
        'INDEX must be a positive integer.');
end
spectral_options.index = index;
[lambda,dlambda,ddlambda_lower,ddlambda_upper,diagnostics] = ...
    calc_ddlambda1_bernstein_strong_bounds( ...
        triangle,e_direction,degree,spectral_options);
end
