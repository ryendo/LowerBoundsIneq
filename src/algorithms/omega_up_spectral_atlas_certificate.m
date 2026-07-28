function certificate = omega_up_spectral_atlas_certificate( ...
    atlas_source,x,y,expected_sha256)
%OMEGA_UP_SPECTRAL_ATLAS_CERTIFICATE  Select one reusable CR/LG anchor.

global INTERVAL_MODE
if isempty(INTERVAL_MODE)
    INTERVAL_MODE = 0;
end

if isstruct(atlas_source)
    atlas = atlas_source;
    atlas_path = '';
    atlas_sha256 = '';
else
    atlas_path = char(atlas_source);
    if nargin < 4 || isempty(expected_sha256)
        error('omega_up_spectral_atlas_certificate:MissingHash', ...
            'A file-backed atlas requires its expected SHA-256 digest.');
    end
    expected_sha256 = lower(char(expected_sha256));
    persistent cached_path cached_sha256 cached_atlas
    if isempty(cached_path) || ~strcmp(cached_path,atlas_path) ...
            || ~strcmp(cached_sha256,expected_sha256)
        atlas_sha256 = ver10_file_sha256(atlas_path);
        if ~strcmp(atlas_sha256,expected_sha256)
            error('omega_up_spectral_atlas_certificate:HashMismatch', ...
                'The spectral-atlas SHA-256 digest does not match.');
        end
        payload = load(atlas_path,'atlas');
        if ~isfield(payload,'atlas')
            error('omega_up_spectral_atlas_certificate:MissingAtlas', ...
                'The atlas MAT file does not contain variable ATLAS.');
        end
        atlas_sha256_after_load = ver10_file_sha256(atlas_path);
        if ~strcmp(atlas_sha256_after_load,expected_sha256)
            error('omega_up_spectral_atlas_certificate:HashMismatch', ...
                'The spectral-atlas file changed while it was loaded.');
        end
        cached_path = atlas_path;
        cached_sha256 = atlas_sha256;
        cached_atlas = payload.atlas;
    else
        atlas_sha256 = cached_sha256;
    end
    atlas = cached_atlas;
end

required = {'schema','rigorous', ...
    'x_domain_lower','x_domain_upper', ...
    'y_domain_lower','y_domain_upper', ...
    'x_anchors','y_anchors','certificates'};
for k = 1:numel(required)
    if ~isfield(atlas,required{k})
        error('omega_up_spectral_atlas_certificate:BadAtlas', ...
            'Missing atlas.%s.',required{k});
    end
end
if ~strcmp(char(atlas.schema), ...
        'lowerboundsineq.omega-up-spectral-cr-lg-atlas.v1')
    error('omega_up_spectral_atlas_certificate:BadAtlas', ...
        'Unsupported spectral-atlas schema.');
end
domain = [atlas.x_domain_lower,atlas.x_domain_upper, ...
    atlas.y_domain_lower,atlas.y_domain_upper];
if ~isnumeric(domain) || numel(domain) ~= 4 ...
        || any(~isfinite(double(domain))) ...
        || double(atlas.x_domain_lower) > double(atlas.x_domain_upper) ...
        || double(atlas.y_domain_lower) > double(atlas.y_domain_upper)
    error('omega_up_spectral_atlas_certificate:BadAtlas', ...
        'The atlas domain endpoints must be finite ordered scalars.');
end
if ~local_boolean_scalar(atlas.rigorous)
    error('omega_up_spectral_atlas_certificate:BadAtlas', ...
        'atlas.rigorous must be scalar 0 or 1.');
end
atlas_rigorous = local_true_scalar(atlas.rigorous);
if INTERVAL_MODE && ~atlas_rigorous
    error('omega_up_spectral_atlas_certificate:UnverifiedAtlas', ...
        'A rigorous Hessian run requires a rigorous spectral atlas.');
end
x_anchors = double(atlas.x_anchors(:).');
y_anchors = double(atlas.y_anchors(:).');
if isempty(x_anchors) || isempty(y_anchors) ...
        || any(~isfinite(x_anchors)) || any(~isfinite(y_anchors)) ...
        || ~isequal(size(atlas.certificates), ...
            [numel(y_anchors),numel(x_anchors)])
    error('omega_up_spectral_atlas_certificate:BadAtlas', ...
        'The atlas anchors and certificate array are incoherent.');
end

x = I_intval(x);
y = I_intval(y);
if I_inf(x) < double(atlas.x_domain_lower) ...
        || I_sup(x) > double(atlas.x_domain_upper) ...
        || I_inf(y) < double(atlas.y_domain_lower) ...
        || I_sup(y) > double(atlas.y_domain_upper)
    error('omega_up_spectral_atlas_certificate:OutsideAtlas', ...
        'The requested cell is outside the certified atlas domain.');
end
[~,ix] = min(abs(x_anchors-I_mid(x)));
[~,iy] = min(abs(y_anchors-I_mid(y)));
if iy > size(atlas.certificates,1) ...
        || ix > size(atlas.certificates,2)
    error('omega_up_spectral_atlas_certificate:BadAtlas', ...
        'The atlas anchor array and certificate array are incoherent.');
end
certificate = atlas.certificates{iy,ix};
local_validate_certificate(certificate);
if ~local_boolean_scalar(certificate.rigorous)
    error('omega_up_spectral_atlas_certificate:BadCertificate', ...
        'certificate.rigorous must be scalar 0 or 1.');
end
if INTERVAL_MODE && ~local_true_scalar(certificate.rigorous)
    error('omega_up_spectral_atlas_certificate:UnverifiedCertificate', ...
        'The selected spectral certificate is not rigorous.');
end
certificate.atlas_index_x = ix;
certificate.atlas_index_y = iy;
certificate.atlas_path = atlas_path;
certificate.atlas_sha256 = atlas_sha256;
end


function local_validate_certificate(certificate)
required = {'schema','reference_triangle','lambda1_LG_lower', ...
    'lambda1_Ritz_upper','lambda2_CR_lower','rigorous'};
for k = 1:numel(required)
    if ~isfield(certificate,required{k})
        error('omega_up_spectral_atlas_certificate:BadCertificate', ...
            'Missing certificate.%s.',required{k});
    end
end
if ~strcmp(char(certificate.schema), ...
        'lowerboundsineq.triangle-spectral-cr-lg-certificate.v1')
    error('omega_up_spectral_atlas_certificate:BadCertificate', ...
        'Unsupported spectral-certificate schema.');
end
triangle = I_intval(certificate.reference_triangle);
if length(triangle) ~= 6
    error('omega_up_spectral_atlas_certificate:BadCertificate', ...
        'The certificate reference triangle must have length six.');
end
prefix = [0,0,1,0];
for k = 1:4
    if I_inf(triangle(k)) ~= prefix(k) ...
            || I_sup(triangle(k)) ~= prefix(k)
        error('omega_up_spectral_atlas_certificate:BadCertificate', ...
            'The certificate reference triangle is not canonical.');
    end
end
if I_inf(triangle(5)) ~= I_sup(triangle(5)) ...
        || I_inf(triangle(6)) ~= I_sup(triangle(6)) ...
        || ~(I_inf(triangle(6)) > 0)
    error('omega_up_spectral_atlas_certificate:BadCertificate', ...
        'The certificate reference must be a point triangle.');
end
end


function tf = local_true_scalar(value)
tf = local_boolean_scalar(value) ...
    && double(value) == 1;
end


function tf = local_boolean_scalar(value)
tf = isscalar(value) ...
    && (islogical(value) || isnumeric(value)) ...
    && isfinite(double(value)) ...
    && any(double(value) == [0,1]);
end
