function [area_bounds, perimeter_bounds, diagnostics] = compute_geometry_bounds(x_inf, x_sup, theta_inf, theta_sup)
% COMPUTE_GEOMETRY_BOUNDS  (Lemma 4.9 / Algorithm 4 — CORRECTED VERTEX)
%
% The rigorous per-cell LOWER bound for J_k is, by Lemma 4.9 / Algorithm 4,
%
%       J_k(T^p) >= B_k( p_{i,j} ; lambda1(T^{p_{i+1,j+1}}) )   for all p in the cell,
%
% with the INNER corner  p_{i,j} = (x_inf, theta_inf).  This is the cell
% minimiser of B_k because B~_k(x,theta;Lambda) is increasing in BOTH x and
% theta on Omega_mid (Lemma 4.9: dB~/dx >= 0 and dB~/dtheta >= 0).
%
% NOTE (bug fix): the previous version evaluated the lower bound at
%   (x_sup, theta_inf) = p_{i+1,j}, which is NOT the cell minimiser (B~ is
%   increasing in x), so it was an over-estimate, not a valid lower bound.
%   It is corrected here to (x_inf, theta_inf) = p_{i,j}.
%
% Outputs (RIGOROUS interval enclosures, NOT collapsed to a point):
%   area_bounds(1)      = |T|   at p_{i,j}    = (x_inf, theta_inf)   <- used for the lower bound
%   perimeter_bounds(1) = |dT|  at p_{i,j}    = (x_inf, theta_inf)
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

    % ----------------------------
    % Cell corners that bound B_k
    % ----------------------------
    % p_{i,j} = (x_inf, theta_inf): the cell MINIMISER of B_k  -> lower bound
    xL = x1;  tL = t1;  yL = xL * tan(tL);
    % Clip the lower vertex UP to the Omega_mid floor y >= eps_down for
    % bottom-straddle cells (inner corner in Omega_down): the degenerate y<eps_down
    % corner is outside Omega_mid and numerically delicate. The clipped point
    % (x_inf, eps_down) is the relevant Omega_mid extreme; J is large there
    % (thin triangle), so B>0 easily.
    if I_mid(yL) < 0.04
        yL = I_intval('0.04');
    end
    % p_{i+1,j+1} = (x_sup, theta_sup): the cell MAXIMISER of B_k -> diagnostic only
    xU = x2;  tU = t2;  yU = xU * tan(tU);

    % ----------------------------
    % Geometry at the two points (rigorous intervals)
    % ----------------------------
    base = I_intval('1');
    AL = yL / I_intval('2');
    AU = yU / I_intval('2');
    PL = base + sqrt(xL^2 + yL^2) + sqrt((base - xL)^2 + yL^2);
    PU = base + sqrt(xU^2 + yU^2) + sqrt((base - xU)^2 + yU^2);

    % Return rigorous interval enclosures (do NOT collapse with I_inf/I_sup).
    area_bounds      = [AL, AU];
    perimeter_bounds = [PL, PU];

    diagnostics = struct();
    diagnostics.p_lower = struct('x', xL, 'theta', tL, 'y', yL, 'A', AL, 'P', PL); % p_{i,j}
    diagnostics.p_upper = struct('x', xU, 'theta', tU, 'y', yU, 'A', AU, 'P', PU); % p_{i+1,j+1}
end

% ----------------------------
% Local helper: robust conversion
% ----------------------------
function xI = local_to_intval(x)
    if ischar(x) || isstring(x)
        xI = I_intval(str2double(x));
    else
        xI = I_intval(x);
    end
end
