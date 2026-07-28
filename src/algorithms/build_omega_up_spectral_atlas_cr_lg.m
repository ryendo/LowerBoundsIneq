function atlas = build_omega_up_spectral_atlas_cr_lg( ...
    eps_up,Nx_anchor,Ny_anchor,N_LG,N_rho,bound_order)
%BUILD_OMEGA_UP_SPECTRAL_ATLAS_CR_LG  Reusable verified endpoint atlas.
%
% Each anchor stores only
%
%   L1_LG <= lambda_1 <= U1_Ritz < L2_CR <= lambda_2.
%
% A fine Omega_up cell reuses the nearest anchor and transports these
% endpoints by the exact two-by-two affine metric factors.  Thus CR/LG is
% not repeated on every Hessian cell.

global INTERVAL_MODE
if isempty(INTERVAL_MODE)
    INTERVAL_MODE = 0;
end
eps_up = I_intval(eps_up);
if ~(I_inf(eps_up) > 0)
    error('build_omega_up_spectral_atlas_cr_lg:BadEpsilon', ...
        'EPS_UP must be strictly positive.');
end
integers = [Nx_anchor,Ny_anchor,N_LG,N_rho,bound_order];
if any(~isfinite(integers)) || any(integers < 1) ...
        || any(integers ~= floor(integers))
    error('build_omega_up_spectral_atlas_cr_lg:BadParameters', ...
        'All mesh and atlas parameters must be positive integers.');
end

x_left = I_intval('0.5');
x_right = x_left+2*eps_up;
y_top = sqrt(I_intval('3'))/2;
y_bottom = y_top-eps_up;
if Nx_anchor == 1
    x_anchors = I_mid((x_left+x_right)/2);
else
    x_anchors = linspace(I_mid(x_left),I_mid(x_right),Nx_anchor);
end
if Ny_anchor == 1
    y_anchors = I_mid((y_bottom+y_top)/2);
else
    y_anchors = linspace(I_mid(y_bottom),I_mid(y_top),Ny_anchor);
end
certificates = cell(Ny_anchor,Nx_anchor);

for iy = 1:Ny_anchor
    for ix = 1:Nx_anchor
        triangle = I_intval( ...
            [0,0,1,0,x_anchors(ix),y_anchors(iy)]);
        certificates{iy,ix} = ...
            triangle_spectral_certificate_cr_lg( ...
                triangle,N_LG,N_rho,bound_order);
        fprintf('spectral atlas: %d/%d anchors complete\n', ...
            (iy-1)*Nx_anchor+ix,Nx_anchor*Ny_anchor);
    end
end

certificate_rigorous = cellfun( ...
    @(c) isfield(c,'rigorous') ...
        && isscalar(c.rigorous) ...
        && (islogical(c.rigorous) || isnumeric(c.rigorous)) ...
        && isfinite(double(c.rigorous)) ...
        && double(c.rigorous) == 1,certificates);
if INTERVAL_MODE && ~all(certificate_rigorous(:))
    error('build_omega_up_spectral_atlas_cr_lg:UnverifiedCertificate', ...
        'At least one anchor is not a rigorous interval certificate.');
end

atlas = struct();
atlas.schema = ...
    'lowerboundsineq.omega-up-spectral-cr-lg-atlas.v1';
atlas.method = ...
    'nearest-anchor-affine-transport-of-LG-Ritz-CR-endpoints';
atlas.rigorous = logical(INTERVAL_MODE && all(certificate_rigorous(:)));
atlas.eps_up = eps_up;
atlas.x_domain = I_infsup(I_inf(x_left),I_sup(x_right));
atlas.y_domain = I_infsup(I_inf(y_bottom),I_sup(y_top));
atlas.x_domain_lower = I_inf(x_left);
atlas.x_domain_upper = I_sup(x_right);
atlas.y_domain_lower = I_inf(y_bottom);
atlas.y_domain_upper = I_sup(y_top);
atlas.x_anchors = x_anchors;
atlas.y_anchors = y_anchors;
atlas.certificates = certificates;
atlas.Nx_anchor = Nx_anchor;
atlas.Ny_anchor = Ny_anchor;
atlas.N_LG = N_LG;
atlas.N_rho = N_rho;
atlas.bound_order = bound_order;
atlas.num_certificates = Nx_anchor*Ny_anchor;
end
