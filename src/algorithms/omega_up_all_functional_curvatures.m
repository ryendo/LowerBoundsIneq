function bounds = omega_up_all_functional_curvatures( ...
    direction,x,y,dlambda,ddlambda_lower,ddlambda_upper)
%OMEGA_UP_ALL_FUNCTIONAL_CURVATURES  Simultaneous J1/J2/Jup certificates.
%
% For T(x,y)=conv{(0,0),(1,0),(x,y)}, set A=y/2 and let P denote the
% perimeter.  The three scale-invariant functionals relevant to Omega_up
% differ only in their explicit geometric term:
%
%   J1  = A*lambda_1-(pi^2/16)P^2/A-constant,
%   J2  = A*lambda_1-C_* (P+sqrt(4*pi*A))^2/(4A),
%   Jup = A*lambda_1-(pi^2/12)P^2/A-constant.
%
% DIRECTION is 'x' or 'y'.  One enclosure of lambda_xx (or lambda_yy),
% together with lambda_y in the y case, therefore certifies all three
% functionals.  J1 and J2 require lower curvature bounds; Jup requires an
% upper curvature bound.

direction = lower(char(direction));
if ~any(strcmp(direction,{'x','y'}))
    error('omega_up_all_functional_curvatures:BadDirection', ...
        'direction must be ''x'' or ''y''.');
end

x = I_intval(x);
y = I_intval(y);
dlambda = I_intval(dlambda);
ddlambda_lower = I_intval(ddlambda_lower);
ddlambda_upper = I_intval(ddlambda_upper);
if ~(I_inf(y) > 0)
    error('omega_up_all_functional_curvatures:BadHeight', ...
        'The triangle height must be certified positive.');
end

if strcmp(direction,'x')
    R1 = local_Rxx('J1',x,y);
    R2 = local_Rxx('J2',x,y);
    Rup = local_Rxx('JUP',x,y);
    expr1 = I_intval('0.5')*y*ddlambda_lower+R1;
    expr2 = I_intval('0.5')*y*ddlambda_lower+R2;
    exprup = I_intval('0.5')*y*ddlambda_upper+Rup;
else
    R1 = local_Ryy('J1',x,y);
    R2 = local_Ryy('J2',x,y);
    Rup = local_Ryy('JUP',x,y);
    expr1 = I_intval('0.5')*y*ddlambda_lower+dlambda+R1;
    expr2 = I_intval('0.5')*y*ddlambda_lower+dlambda+R2;
    exprup = I_intval('0.5')*y*ddlambda_upper+dlambda+Rup;
end

bounds = struct();
bounds.direction = direction;
bounds.J1 = local_lower_certificate(expr1,R1);
bounds.J2 = local_lower_certificate(expr2,R2);
bounds.JUP = local_upper_certificate(exprup,Rup);
bounds.all_signs_verified = ...
    bounds.J1.sign_verified ...
    && bounds.J2.sign_verified ...
    && bounds.JUP.sign_verified;
end


function certificate = local_lower_certificate(expression,remainder)
certificate = struct();
certificate.sense = 'lower-positive';
certificate.expression = expression;
certificate.remainder = remainder;
certificate.bound = I_intval(I_inf(expression));
bound_lower = I_inf(certificate.bound);
certificate.sign_verified = isfinite(bound_lower) && bound_lower > 0;
end


function certificate = local_upper_certificate(expression,remainder)
certificate = struct();
certificate.sense = 'upper-negative';
certificate.expression = expression;
certificate.remainder = remainder;
certificate.bound = I_intval(I_sup(expression));
bound_upper = I_sup(certificate.bound);
certificate.sign_verified = isfinite(bound_upper) && bound_upper < 0;
end


function Rxx = local_Rxx(kind,x,y)
r1 = sqrt(x.^2+y.^2);
r2 = sqrt((x-I_intval(1)).^2+y.^2);
P = I_intval(1)+r1+r2;
geometry = (x./r1+(x-I_intval(1))./r2).^2 ...
    +P.*y.^2.*(1./r1.^3+1./r2.^3);

switch kind
    case 'J1'
        Rxx = -(I_pi^2./(4*y)).*geometry;
    case 'JUP'
        Rxx = -(I_pi^2./(3*y)).*geometry;
    case 'J2'
        Cstar = 4*I_pi^2 ...
            /(3+sqrt(I_pi*sqrt(I_intval(3))))^2;
        Q = P+sqrt(2*I_pi*y);
        geometry2 = (x./r1+(x-I_intval(1))./r2).^2 ...
            +Q.*y.^2.*(1./r1.^3+1./r2.^3);
        Rxx = -(Cstar./y).*geometry2;
    otherwise
        error('omega_up_all_functional_curvatures:InternalFunctional', ...
            'Unknown functional %s.',kind);
end
end


function Ryy = local_Ryy(kind,x,y)
r1 = sqrt(x.^2+y.^2);
r2 = sqrt((x-I_intval(1)).^2+y.^2);
P = I_intval(1)+r1+r2;
geometry = y.^2.*P.*(x.^2./r1.^3 ...
    +(x-I_intval(1)).^2./r2.^3) ...
    +(I_intval(1)+x.^2./r1 ...
      +(x-I_intval(1)).^2./r2).^2;

switch kind
    case 'J1'
        Ryy = -(I_pi^2./(4*y.^3)).*geometry;
    case 'JUP'
        Ryy = -(I_pi^2./(3*y.^3)).*geometry;
    case 'J2'
        Cstar = 4*I_pi^2 ...
            /(3+sqrt(I_pi*sqrt(I_intval(3))))^2;
        Q = P+sqrt(2*I_pi*y);
        Qy = y./r1+y./r2+sqrt(I_pi./(2*y));
        Qyy = x.^2./r1.^3 ...
            +(x-I_intval(1)).^2./r2.^3 ...
            -I_intval('0.5')*sqrt(I_pi/2).*y.^(-3/2);
        Ryy = -(Cstar/2).*( ...
            (2./y).*(Qy.^2+Q.*Qyy) ...
            -(4./y.^2).*(Q.*Qy) ...
            +(2./y.^3).*Q.^2);
    otherwise
        error('omega_up_all_functional_curvatures:InternalFunctional', ...
            'Unknown functional %s.',kind);
end
end
