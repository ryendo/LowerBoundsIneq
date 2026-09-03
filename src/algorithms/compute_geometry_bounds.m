function [area_bounds, perimeter_bounds, diagnostics] = compute_geometry_bounds(x_inf, x_sup, theta_inf, theta_sup)
% COMPUTE_GEOMETRY_BOUNDS  (Lemma 4.9 / Algorithm 4 — COMPARISON POINT)
%
% The rigorous per-cell lower bound has the form
%
%   J_k(T^p) >= B_k(q_C;lambda1(T^{p_{i+1,j+1}})).
%
% If the lower-left parameter corner is certified to lie at or above
% y=0.04, q_C is that corner.  Otherwise q_C=(x_sup,0.04).  The latter is
% an auxiliary comparison point, not necessarily a minimizer in the target
% portion and not necessarily a point of Omega (x_sup=1 is allowed).  For
% every target p=(x,y) and Lambda>=0, monotonicity in y at fixed x and
% monotonicity of the perimeter at fixed y give
%
%   B_k((x_sup,0.04);Lambda) <= B_k((x,0.04);Lambda)
%                              <= B_k(p;Lambda).
%
% Outputs (RIGOROUS interval enclosures, NOT collapsed to a point):
%   area_bounds(1), perimeter_bounds(1) are evaluated at the certified
%   comparison point.  Usually this is p_{i,j}=(x_inf,theta_inf).  If the
%   lower-left corner is not certified above y=0.04, the comparison point
%   is instead (x_sup,0.04); see the proof above.
%   area_bounds(2)      = |T|   at p_{i+1,j+1}= (x_sup, theta_sup)   <- diagnostic (upper vertex)
%   perimeter_bounds(2) = |dT|  at p_{i+1,j+1}= (x_sup, theta_sup)
%
% The interval enclosures are propagated all the way into B_k so that
% I_inf(B_k) is a rigorous lower bound (the geometry intervals are thin but
% B_k is not monotone in A and P separately, so we must NOT collapse them).

    % ----------------------------
    % Convert inputs to interval
    % ----------------------------
    x1 = local_to_intval(x_inf);
    x2 = local_to_intval(x_sup);
    t1 = local_to_intval(theta_inf);
    t2 = local_to_intval(theta_sup);
    local_require_finite_scalar_interval(x1,'x_inf');
    local_require_finite_scalar_interval(x2,'x_sup');
    local_require_finite_scalar_interval(t1,'theta_inf');
    local_require_finite_scalar_interval(t2,'theta_sup');

    % ----------------------------
    % Cell corners that bound B_k
    % ----------------------------
    % Candidate lower-left comparison point p_{i,j}=(x_inf,theta_inf).
    xCorner = x1;
    tCorner = t1;
    yCorner = xCorner * tan(tCorner);
    yFloor = I_intval('0.04');

    % On Omega, B_k(x,theta;Lambda) is nondecreasing in each of x and
    % theta.  When the lower-left corner is not certified above the floor,
    % simply clipping it to (x_inf,0.04) is not a lower bound: along the
    % horizontal line y=0.04, the perimeter (and hence the subtracted
    % geometric term) increases with x.  For every target point p in the
    % cell with y>=0.04,
    %
    %   B_k(p;Lambda) >= B_k(x_p,0.04;Lambda)
    %                 >= B_k(x_sup,0.04;Lambda).
    %
    % Use this right-hand floor point unless the lower-left corner is
    % rigorously certified above the floor.  An interval overlap therefore
    % takes the conservative floor branch.  Reject nonfinite arithmetic
    % explicitly: NaN comparisons in MATLAB are false and hence are not by
    % themselves fail-closed.
    local_require_finite_scalar_interval( ...
        yCorner,'lower-left y coordinate');
    certifiedAbove = I_inf(yCorner) >= I_sup(yFloor);
    useFloorComparison = ~certifiedAbove;
    if useFloorComparison
        xL = x2;
        tL = atan(yFloor/xL);
        yL = yFloor;
        comparisonRule = 'right endpoint on y floor comparison point';
    else
        xL = xCorner;
        tL = tCorner;
        yL = yCorner;
        comparisonRule = 'lower-left parameter comparison point';
    end
    % p_{i+1,j+1}=(x_sup,theta_sup): outer spectral corner and diagnostic.
    xU = x2;  tU = t2;  yU = xU * tan(tU);
    local_require_finite_scalar_interval(yU,'upper-right y coordinate');
    cellCertifiedBelowFloor = I_sup(yU) < I_inf(yFloor);

    % ----------------------------
    % Geometry at the two points (rigorous intervals)
    % ----------------------------
    base = I_intval('1');
    AL = yL / I_intval('2');
    AU = yU / I_intval('2');
    PL = base + sqrt(xL^2 + yL^2) + sqrt((base - xL)^2 + yL^2);
    PU = base + sqrt(xU^2 + yU^2) + sqrt((base - xU)^2 + yU^2);
    local_require_finite_scalar_interval(AL,'comparison-point area');
    local_require_finite_scalar_interval( ...
        PL,'comparison-point perimeter');
    local_require_finite_scalar_interval(AU,'upper-right area');
    local_require_finite_scalar_interval(PU,'upper-right perimeter');

    % Return rigorous interval enclosures (do NOT collapse with I_inf/I_sup).
    area_bounds      = [AL, AU];
    perimeter_bounds = [PL, PU];

    diagnostics = struct();
    diagnostics.p_comparison = struct( ...
        'x',xL,'theta',tL,'y',yL,'A',AL,'P',PL);
    % Backward-compatible alias; this is a comparison point, not a claimed
    % minimizer of the target portion.
    diagnostics.p_lower = diagnostics.p_comparison;
    diagnostics.p_upper = struct('x', xU, 'theta', tU, 'y', yU, 'A', AU, 'P', PU); % p_{i+1,j+1}
    diagnostics.lower_corner_certified_above_floor = certifiedAbove;
    diagnostics.use_floor_comparison = useFloorComparison;
    diagnostics.cell_certified_below_floor = cellCertifiedBelowFloor;
    diagnostics.comparison_geometry_rule = comparisonRule;
    % Backward-compatible aliases for existing diagnostic consumers.
    diagnostics.bottom_straddle = useFloorComparison;
    diagnostics.lower_geometry_rule = comparisonRule;
end

% ----------------------------
% Local helper: robust conversion
% ----------------------------
function xI = local_to_intval(x)
    if ischar(x) || isstring(x)
        % Preserve the exact decimal endpoint from cell_def.csv.  Passing
        % through a binary double first can move a one-sided cell boundary.
        xI = I_intval(char(x));
    else
        xI = I_intval(x);
    end
end


function local_require_finite_scalar_interval(value,name)
    lower = I_inf(value);
    upper = I_sup(value);
    if ~isscalar(lower) || ~isscalar(upper) ...
            || ~isfinite(lower) || ~isfinite(upper)
        error('compute_geometry_bounds:NonfiniteGeometry', ...
            '%s must have finite scalar interval endpoints.',name);
    end
end
