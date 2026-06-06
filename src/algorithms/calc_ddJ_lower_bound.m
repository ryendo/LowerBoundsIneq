function [ddJ_lower, components] = calc_ddJ_lower_bound(conjecture_type, triangle, e_direction, N_spectral, mesh_params)
    % CALC_DDJ_LOWER_BOUND: Compute lower bound for second-order directional derivative of J1 or J2
    %
    % Paper Reference: Section 3 (Main Theory) and Section 4 (Computer-assisted proof)
    %   The second-order directional derivative of J = f(lambda_1, Area, Perimeter) is
    %   computed using the chain rule and rigorous bounds for each component.
    %
    % For J1(triangle) = lambda_1 * Area - C1 * Perim^2 / Area - C2:
    %   ddJ1/de^2 = ddlam1 * Area + 2*dlam1 * dArea + lam1 * ddArea
    %             - C1 * d^2/de^2[Perim^2/Area]
    %
    % For the equilateral triangle with direction e = (a, b):
    %   - dArea/de = b/2 (Area = y/2 = x*tan(theta)/2, and y depends on b)
    %   - ddArea/de^2 = 0 (Area is linear in y for fixed base)
    %   - dPerim/de, ddPerim/de^2 depend on side lengths
    %
    % Inputs:
    %   conjecture_type: 'J1' (Laugesen-Siudeja) or 'J2' (Cheeger-type)
    %   triangle: [x1 y1 x2 y2 x3 y3] triangle vertices (row vector)
    %   e_direction: [a, b] direction vector for directional derivative
    %   N_spectral: number of terms in spectral truncation (for ddlam1)
    %   mesh_params: struct with fields N_LG, N_rho, fem_ord, fem_ord_LG
    %
    % Outputs:
    %   ddJ_lower: rigorous lower bound on second-order directional derivative of J
    %   components: detailed breakdown of computation
    %
    % Author: Based on paper by R. Endo, X. Liu, P. Mariano
    % Date: 2025-01-14

    fprintf('=== calc_ddJ_lower_bound for %s ===\n', conjecture_type);

    %% Extract triangle geometry
    % Triangle vertices: (0,0), (1,0), (x0, y0)
    x0 = triangle(5);
    y0 = triangle(6);

    % Direction vector
    a = e_direction(1);
    b = e_direction(2);

    % Convert to interval
    x0_intval = I_intval(x0);
    y0_intval = I_intval(y0);
    a_intval = I_intval(a);
    b_intval = I_intval(b);

    %% Step 1: Compute second-order derivative of lambda_1
    fprintf('Step 1: Computing ddlam1 lower bound...\n');

    % Create base triangle for calc_ddlami_lower_bound
    base_triangle = I_intval([0, 0, 1, 0, x0, y0]);  % Same as target for equilateral case

    % Call existing function with i=1 (lambda_1)
    [lam1, dlam1, ddlam1_lower] = calc_ddlami_lower_bound(1, base_triangle, triangle, ...
        e_direction, N_spectral, mesh_params.N_LG, mesh_params.N_rho, ...
        mesh_params.fem_ord, mesh_params.fem_ord_LG);

    fprintf('  lam1 = %.17f\n', I_mid(lam1));
    fprintf('  dlam1 = %.17f\n', I_mid(dlam1));
    fprintf('  ddlam1_lower = %.17f\n', I_inf(ddlam1_lower));

    %% Step 2: Compute geometry and its derivatives
    fprintf('Step 2: Computing geometry derivatives...\n');

    % Current geometry (at p0)
    % Area = y0/2 (for triangle with base 1 and apex at (x0, y0))
    Area = y0_intval / I_intval('2');

    % Perimeter = 1 + sqrt(x0^2 + y0^2) + sqrt((1-x0)^2 + y0^2)
    side1 = I_intval('1');
    side2 = sqrt(x0_intval^2 + y0_intval^2);
    side3 = sqrt((I_intval('1') - x0_intval)^2 + y0_intval^2);
    Perim = side1 + side2 + side3;

    %% Step 2.1: First-order derivatives of geometry
    % For perturbation p -> p + t*e = (x0 + t*a, y0 + t*b):
    %
    % dArea/dt = d(y/2)/dt = b/2
    dArea = b_intval / I_intval('2');

    % dPerim/dt = d(side2)/dt + d(side3)/dt
    %   side2 = sqrt(x^2 + y^2)
    %   d(side2)/dt = (x*a + y*b) / sqrt(x^2 + y^2) = (x0*a + y0*b) / side2
    d_side2 = (x0_intval * a_intval + y0_intval * b_intval) / side2;

    %   side3 = sqrt((1-x)^2 + y^2)
    %   d(side3)/dt = (-(1-x)*a + y*b) / sqrt((1-x)^2 + y^2) = (-(1-x0)*a + y0*b) / side3
    d_side3 = (-(I_intval('1') - x0_intval) * a_intval + y0_intval * b_intval) / side3;

    dPerim = d_side2 + d_side3;

    fprintf('  dArea = %.10f\n', I_mid(dArea));
    fprintf('  dPerim = %.10f\n', I_mid(dPerim));

    %% Step 2.2: Second-order derivatives of geometry
    % ddArea/dt^2 = 0 (Area is linear in y)
    ddArea = I_intval('0');

    % dd(side2)/dt^2 = d/dt[(x*a + y*b)/sqrt(x^2 + y^2)]
    %   = (a^2 + b^2)/sqrt(x^2+y^2) - (x*a + y*b)^2 / (x^2+y^2)^(3/2)
    %   = [(a^2 + b^2) * (x^2 + y^2) - (x*a + y*b)^2] / (x^2 + y^2)^(3/2)
    %   = [a^2*y^2 + b^2*x^2 - 2*a*b*x*y] / (x^2 + y^2)^(3/2)  [by expansion]
    %   = (a*y - b*x)^2 / (x^2 + y^2)^(3/2)

    numer2 = (a_intval * y0_intval - b_intval * x0_intval)^2;
    denom2 = (x0_intval^2 + y0_intval^2)^(I_intval('3')/I_intval('2'));
    dd_side2 = numer2 / denom2;

    % dd(side3)/dt^2 = (a*y - b*(1-x))^2 / ((1-x)^2 + y^2)^(3/2)
    numer3 = (a_intval * y0_intval - b_intval * (I_intval('1') - x0_intval))^2;
    denom3 = ((I_intval('1') - x0_intval)^2 + y0_intval^2)^(I_intval('3')/I_intval('2'));
    dd_side3 = numer3 / denom3;

    ddPerim = dd_side2 + dd_side3;

    fprintf('  ddArea = %.10f\n', I_mid(ddArea));
    fprintf('  ddPerim = %.10f\n', I_mid(ddPerim));

    %% Step 3: Compute second-order derivative of J
    fprintf('Step 3: Computing ddJ lower bound...\n');

    % Constants
    % Route through I_pi so that, in interval mode, pi is enclosed via
    % intval('pi') (rigorous enclosure of the true constant). I_intval(pi)
    % would enclose only the double-precision value of pi, which does not
    % contain the exact pi.
    pi_val = I_pi;

    if strcmpi(conjecture_type, 'J1')
        % J1 = lambda_1 * Area - C1 * Perim^2 / Area - C2
        % C1 = pi^2/16, C2 = 7*sqrt(3)*pi^2/12
        C1 = pi_val^2 / I_intval('16');

        % ddJ1 = ddlam1*Area + 2*dlam1*dArea + lam1*ddArea
        %      - C1 * dd(Perim^2/Area)
        %
        % Let f = Perim^2/Area
        % df = (2*Perim*dPerim*Area - Perim^2*dArea) / Area^2
        %    = (2*Perim*dPerim)/Area - Perim^2*dArea/Area^2
        % ddf = d/dt[(2*Perim*dPerim)/Area - Perim^2*dArea/Area^2]
        %     = (2*dPerim^2 + 2*Perim*ddPerim)/Area - 2*Perim*dPerim*dArea/Area^2
        %       - (2*Perim*dPerim*dArea)/Area^2 - Perim^2*ddArea/Area^2 + 2*Perim^2*dArea^2/Area^3

        % Simplify ddf:
        % ddf = 2*dPerim^2/Area + 2*Perim*ddPerim/Area - 4*Perim*dPerim*dArea/Area^2
        %       - Perim^2*ddArea/Area^2 + 2*Perim^2*dArea^2/Area^3

        term_ddf_1 = I_intval('2') * dPerim^2 / Area;
        term_ddf_2 = I_intval('2') * Perim * ddPerim / Area;
        term_ddf_3 = -I_intval('4') * Perim * dPerim * dArea / Area^2;
        term_ddf_4 = -Perim^2 * ddArea / Area^2;
        term_ddf_5 = I_intval('2') * Perim^2 * dArea^2 / Area^3;

        ddf = term_ddf_1 + term_ddf_2 + term_ddf_3 + term_ddf_4 + term_ddf_5;

        % ddJ1 = ddlam1*Area + 2*dlam1*dArea + lam1*ddArea - C1*ddf
        term_J_1 = ddlam1_lower * Area;  % Use lower bound
        term_J_2 = I_intval('2') * dlam1 * dArea;
        term_J_3 = lam1 * ddArea;
        term_J_4 = -C1 * ddf;

        ddJ_interval = term_J_1 + term_J_2 + term_J_3 + term_J_4;
        ddJ_lower = I_inf(ddJ_interval);

        components.conjecture = 'J1';
        components.C1 = C1;
        components.ddf = ddf;
        components.term_J = [term_J_1, term_J_2, term_J_3, term_J_4];

    elseif strcmpi(conjecture_type, 'J2')
        % J2 = lambda_1 * Area - C_cheeger * (Perim + sqrt(4*pi*Area))^2 / (4*Area)
        % This is more complex due to the sqrt(Area) term

        sqrt3 = sqrt(I_intval('3'));
        inner_sqrt = sqrt(pi_val * sqrt3);
        denom_c = (I_intval('3') + inner_sqrt)^2;
        C_cheeger = I_intval('4') * pi_val^2 / denom_c;

        % Let g = (Perim + sqrt(4*pi*Area))^2 / (4*Area)
        % Let h = sqrt(4*pi*Area) = sqrt(4*pi) * sqrt(Area) = 2*sqrt(pi*Area)
        % dh = sqrt(pi/Area) * dArea = sqrt(pi) * dArea / sqrt(Area)
        % ddh = sqrt(pi) * [ddArea/sqrt(Area) - dArea^2/(2*Area^(3/2))]
        %     = sqrt(pi) * [ddArea*sqrt(Area) - dArea^2/(2*sqrt(Area))] / Area

        sqrt_pi = sqrt(pi_val);
        sqrt_Area = sqrt(Area);
        h = I_intval('2') * sqrt(pi_val * Area);
        dh = sqrt_pi * dArea / sqrt_Area;
        ddh = sqrt_pi * (ddArea / sqrt_Area - dArea^2 / (I_intval('2') * Area * sqrt_Area));

        % Numerator: N = (Perim + h)^2
        % dN = 2*(Perim + h)*(dPerim + dh)
        % ddN = 2*(dPerim + dh)^2 + 2*(Perim + h)*(ddPerim + ddh)
        N = (Perim + h)^2;
        dN = I_intval('2') * (Perim + h) * (dPerim + dh);
        ddN = I_intval('2') * (dPerim + dh)^2 + I_intval('2') * (Perim + h) * (ddPerim + ddh);

        % Denominator: D = 4*Area
        % dD = 4*dArea
        % ddD = 4*ddArea = 0
        D = I_intval('4') * Area;
        dD = I_intval('4') * dArea;
        ddD = I_intval('4') * ddArea;

        % g = N/D
        % dg = (dN*D - N*dD) / D^2
        % ddg = (ddN*D^2 - 2*dN*D*dD - N*D*ddD + 2*N*dD^2) / D^3
        %     = (ddN*D - 2*dN*dD - N*ddD + 2*N*dD^2/D) / D^2

        ddg = (ddN * D - I_intval('2') * dN * dD - N * ddD + I_intval('2') * N * dD^2 / D) / D^2;

        % ddJ2 = ddlam1*Area + 2*dlam1*dArea + lam1*ddArea - C_cheeger*ddg
        term_J_1 = ddlam1_lower * Area;
        term_J_2 = I_intval('2') * dlam1 * dArea;
        term_J_3 = lam1 * ddArea;
        term_J_4 = -C_cheeger * ddg;

        ddJ_interval = term_J_1 + term_J_2 + term_J_3 + term_J_4;
        ddJ_lower = I_inf(ddJ_interval);

        components.conjecture = 'J2';
        components.C_cheeger = C_cheeger;
        components.ddg = ddg;
        components.term_J = [term_J_1, term_J_2, term_J_3, term_J_4];

    else
        error('Unknown conjecture type: %s. Use ''J1'' or ''J2''.', conjecture_type);
    end

    fprintf('  ddJ_lower = %.17f\n', ddJ_lower);

    %% Store all components for diagnostics
    components.triangle = triangle;
    components.e_direction = e_direction;
    components.x0 = x0;
    components.y0 = y0;
    components.Area = Area;
    components.Perim = Perim;
    components.dArea = dArea;
    components.dPerim = dPerim;
    components.ddArea = ddArea;
    components.ddPerim = ddPerim;
    components.lam1 = lam1;
    components.dlam1 = dlam1;
    components.ddlam1_lower = ddlam1_lower;
    components.ddJ_interval = ddJ_interval;
    components.ddJ_lower = ddJ_lower;

    fprintf('=== calc_ddJ_lower_bound complete ===\n');

end
