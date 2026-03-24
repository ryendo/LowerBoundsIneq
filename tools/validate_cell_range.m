% Strict coverage test of Omega_mid by cells in ../inputs/cell_def_v4.csv
clear; clc;

% Parameters of Omega_mid
epsUp   = I_intval('0.122');
epsDown = I_intval('0.04');
xMin    = I_intval('0.5');
xMax    = I_intval('1.0');
tol     = I_intval('1.0e-12');

% Read cells defined in (x, theta)
T  = readtable('inputs/cell_def.csv');
xL = T.x_inf;
xR = T.x_sup;
tL = T.theta_inf;
tR = T.theta_sup;

% Convert interval constants to numeric endpoints where MATLAB needs doubles
xMin_d = inf(xMin);
xMax_d = sup(xMax);
tol_d  = sup(tol);

% Critical x where the top boundary of Omega meets the lower branch
% of the excluded ball around (1/2, sqrt(3)/2)
a = I_intval('2') - epsUp^2;
xSwitch = (a + sqrt(max(I_intval('0'), I_intval('3') * (I_intval('4') - a^2)))) / I_intval('4');

% Critical x where theta_ball(x) attains its minimum
dCrit = I_intval('0.5') * (-epsUp^2 + epsUp * sqrt(max(I_intval('0'), I_intval('3') * (I_intval('1') - epsUp^2))));
xBallMin = I_intval('0.5') + dCrit;

% Build x-slabs:
% inside each slab, the active set of cells is constant,
% and the upper theta-bound of Omega_mid is monotone
bp = [ ...
    xMin_d; ...
    xMax_d; ...
    xL; ...
    xR; ...
    inf(xSwitch); sup(xSwitch); ...
    inf(xBallMin); sup(xBallMin); ...
    inf(I_intval('0.5') + epsUp); sup(I_intval('0.5') + epsUp) ...
];
bp = bp(bp >= xMin_d - tol_d & bp <= xMax_d + tol_d);
bp = sort(bp);
bp = bp(bp >= xMin_d - tol_d & bp <= xMax_d + tol_d);
bp = bp([true; diff(bp) > tol_d]);

isCovered = true;
witness = [];

for k = 1:numel(bp) - 1
    xa = bp(k);
    xb = bp(k+1);

    if xb - xa <= tol_d
        continue;
    end

    % Cells that cover the whole x-slab [xa, xb]
    active = (xL <= xa + tol_d) & (xR >= xb - tol_d);

    % Required theta interval for Omega_mid on this slab:
    % theta >= atan(epsDown / x)
    % theta <= thetaMidUpper(x)
    %
    % Since thetaLow(x) is decreasing, its worst case is at x = xb.
    % Since thetaUp(x) is monotone on each slab, its worst case is at an endpoint.
    thNeedL = inf(thetaMidLower(xb, epsDown));
    thNeedU = sup(max(thetaMidUpper(xa, epsUp), thetaMidUpper(xb, epsUp)));

    % Empty target slab => nothing to check
    if thNeedU < thNeedL - tol_d
        continue;
    end

    % No active cell but target is nonempty => uncovered
    if ~any(active)
        isCovered = false;
        witness = [xa, xb, thNeedL, thNeedU];
        break;
    end

    % Merge theta-intervals of active cells
    intervals = sortrows([tL(active), tR(active)], 1);
    merged = mergeIntervals(intervals, tol_d);

    % Check whether one merged interval contains the whole required theta-range
    ok = any(merged(:,1) <= thNeedL + tol_d & merged(:,2) >= thNeedU - tol_d);

    if ~ok
        isCovered = false;
        witness = [xa, xb, thNeedL, thNeedU];
        break;
    end
end

if isCovered
    fprintf('OK: ../inputs/cell_def.csv completely covers Omega_mid.\n');
else
    fprintf('NG: coverage fails on x in [%.16f, %.16f].\n', witness(1), witness(2));
    fprintf('Required theta interval on this slab is [%.16f, %.16f].\n', witness(3), witness(4));
end

% ===== local functions =====

function th = thetaMidLower(x, epsDown)
    % Lower boundary of Omega_mid: y >= epsDown
    % In (x, theta), this is theta >= atan(epsDown / x)
    th = atan(epsDown ./ x);
end

function th = thetaMidUpper(x, epsUp)
    % Upper boundary of Omega_mid in (x, theta)
    %
    % Omega itself gives y <= sqrt(1 - x^2).
    % Outside the excluded ball around (1/2, sqrt(3)/2) means
    % y <= sqrt(3)/2 - sqrt(epsUp^2 - (x - 1/2)^2)
    % whenever the square root is real.
    %
    % Therefore the upper y-bound is the minimum of these two curves.
    yOmega = sqrt(max(I_intval('0'), I_intval('1') - x.^2));

    r2 = epsUp^2 - (x - I_intval('0.5')).^2;
    yBallLower = I_intval('inf');
    if r2 >= I_intval('0')
        yBallLower = sqrt(I_intval('3'))/I_intval('2') - sqrt(r2);
    end

    yUpper = min(yOmega, yBallLower);
    th = atan(yUpper ./ x);
end

function merged = mergeIntervals(intervals, tol)
    % Merge closed intervals [a,b] with tolerance tol
    merged = [];
    if isempty(intervals)
        return;
    end

    cur = intervals(1,:);
    for i = 2:size(intervals,1)
        nxt = intervals(i,:);
        if nxt(1) <= cur(2) + tol
            cur(2) = max(cur(2), nxt(2));
        else
            merged = [merged; cur]; %#ok<AGROW>
            cur = nxt;
        end
    end
    merged = [merged; cur];
end