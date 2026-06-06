function [area_bounds, perimeter_bounds, diagnostics] = compute_geometry_bounds(x_inf, x_sup, theta_inf, theta_sup)
% COMPUTE_GEOMETRY_BOUNDS (Lemma-aligned):
% Return geometry evaluated at the two Lemma vertices:
%   - Lower-bound geometry point: p_{i+1,j}   = (x_sup, theta_inf)
%   - Upper-bound geometry point: p_{i,j+1}   = (x_inf, theta_sup)
%
% Outputs:
%   area_bounds     = [A_at_p_i1j, A_at_p_ij1] (NOT an enclosure interval!)
%   perimeter_bounds= [P_at_p_i1j, P_at_p_ij1] (NOT an enclosure interval!)
%   diagnostics: struct with the two points in (x,y) and (x,theta)

    % ----------------------------
    % Convert inputs to interval
    % ----------------------------
    x1 = local_to_intval(x_inf);
    x2 = local_to_intval(x_sup);
    t1 = local_to_intval(theta_inf);
    t2 = local_to_intval(theta_sup);

    % ----------------------------
    % Lemma vertices
    % ----------------------------
    % p_{i+1,j}: (x_sup, theta_inf)
    xL = x2;
    tL = t1;
    yL = xL * tan(tL);

    % p_{i,j+1}: (x_inf, theta_sup)
    xU = x1;
    tU = t2;
    yU = xU * tan(tU);

    % ----------------------------
    % Geometry at the two points
    % ----------------------------
    AL = yL / I_intval('2');
    AU = yU / I_intval('2');

    base = I_intval('1');
    PL = base + sqrt(xL^2 + yL^2) + sqrt((I_intval('1') - xL)^2 + yL^2);
    PU = base + sqrt(xU^2 + yU^2) + sqrt((I_intval('1') - xU)^2 + yU^2);

    % Return as numeric scalars (inf=sup at a point, but keep rigorous extraction)
    area_bounds      = [I_inf(AL), I_inf(AU)];
    perimeter_bounds = [I_inf(PL), I_inf(PU)];

    diagnostics = struct();
    diagnostics.p_lower = struct('x', xL, 'theta', tL, 'y', yL, 'A', AL, 'P', PL); % p_{i+1,j}
    diagnostics.p_upper = struct('x', xU, 'theta', tU, 'y', yU, 'A', AU, 'P', PU); % p_{i,j+1}
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
