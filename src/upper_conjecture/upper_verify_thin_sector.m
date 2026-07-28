function result = upper_verify_thin_sector(mode, y_max)
%UPPER_VERIFY_THIN_SECTOR Verify the upper conjecture for 0<y<=y_max.
%
% For T(x,y)=conv{(0,0),(1,0),(x,y)} in diameter normalization, the
% Freitas--Siudeja sector inclusion and the sector trial
%
%   r^m(1-r^m) sin(nu*theta),  nu=pi/asin(y), m=(4/5)nu^(2/3),
%
% give
%
%   lambda_1(T) <= 2 Q_nu/(1+sqrt(1-y^2)),
%   Q_nu=(m+1)(2m+1)(3m+2)(nu^2+2m^2)/(6m^3).
%
% Multiplication by y^2 and z=nu^(-1/3) remove the apparent singularity:
%
% y^2 lambda_1(T) <=
%   2*pi^2*cosc(s)*P(z^2/c)*(1+2*c^2*z^2),
%
% where s=asin(y), c=4/5, cosc(s)=(1-cos(s))/s^2 and
% P(w)=(1+w)(2+w)(3+2w)/6.  On 0<=s<=asin(0.06), the alternating
% Taylor enclosure
%
%   1/2-s^2/24 <= cosc(s) <= 1/2-s^2/24+s^4/720
%
% is rigorous.  Since L>=2, it is enough to compare with
% 4*pi^2/3+2*sqrt(3)*pi^2*y/3.

if nargin < 1 || isempty(mode)
    mode = 'double';
end
if nargin < 2 || isempty(y_max)
    y_max = 0.06;
end
if y_max <= 0 || y_max > 0.06
    error('This certificate is implemented only for 0 < y_max <= 0.06.');
end

is_interval = strcmpi(mode, 'interval');
if is_interval
    y = infsup(0, y_max);
    one = intval(1);
    pi_value = intval('pi');
    c = intval(4) / 5;
else
    % In double mode the endpoint interval is emulated by the monotone
    % endpoint choices below.  This is a diagnostic, not a certificate.
    y = [0, y_max];
    one = 1;
    pi_value = pi;
    c = 4/5;
end

if is_interval
    s = asin(y);
    z = (s / pi_value)^(intval(1)/3);
    w = z^2 / c;
    P = (one+w)*(2*one+w)*(3*one+2*w)/6;
    cosc = intval(1)/2 - s^2/24 + infsup(0,1)*s^4/720;
    H_trial = 2*pi_value^2*cosc*P*(one+2*c^2*z^2);
    H_target = 4*pi_value^2/3 + 2*sqrt(intval(3))*pi_value^2*y/3;
    margin = H_target - H_trial;
    [margin_lo, margin_hi] = upper_interval_endpoints(margin, mode);
    [trial_lo, trial_hi] = upper_interval_endpoints(H_trial, mode);
    [target_lo, target_hi] = upper_interval_endpoints(H_target, mode);
else
    % Report the same natural-enclosure diagnostic using ordinary doubles.
    s_hi = asin(y_max);
    z_hi = (s_hi/pi)^(1/3);
    w_hi = z_hi^2/c;
    P_hi = (1+w_hi)*(2+w_hi)*(3+2*w_hi)/6;
    cosc_hi = 1/2 + s_hi^4/720; % safe diagnostic overestimate
    trial_hi = 2*pi^2*cosc_hi*P_hi*(1+2*c^2*z_hi^2);
    trial_lo = pi^2;
    target_lo = 4*pi^2/3;
    target_hi = 4*pi^2/3 + 2*sqrt(3)*pi^2*y_max/3;
    margin_lo = target_lo - trial_hi;
    margin_hi = target_hi - trial_lo;
end

result.region = 'thin_sector';
result.mode = lower(mode);
result.rigor = ternary(is_interval, 'certified_interval', 'exploratory_double');
result.y_min = 0;
result.y_max = y_max;
result.scaled_trial_lower = trial_lo;
result.scaled_trial_upper = trial_hi;
result.scaled_target_lower = target_lo;
result.scaled_target_upper = target_hi;
result.scaled_margin_lower = margin_lo;
result.scaled_margin_upper = margin_hi;
result.verified = is_interval && margin_lo > 0;
end

function value = ternary(condition, yes_value, no_value)
if condition
    value = yes_value;
else
    value = no_value;
end
end
