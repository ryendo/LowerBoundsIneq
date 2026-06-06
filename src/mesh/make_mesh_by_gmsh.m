function mesh = make_mesh_by_gmsh(a, b, h, block_size, verbose, tag)
% MAKE_MESH_BY_GMSH  Generate a mesh on triangle (0,0)-(1,0)-(a,b).
%
% Determinism guarantee:
%   Per-cell reference mesh. The reference geometry passed to
%   make_mesh_by_gmsh_ref is (I_mid(a), I_mid(b)) itself. The persistent
%   cache is indexed by (tag, h_mid, a_mid, b_mid) with full precision, so
%   a cache hit implies bit-identical reference geometry. Cache state for
%   any key depends only on that key (pure memoization) — no call-order
%   dependence. Hence resume and single-pass runs produce identical
%   meshes, and J_lower is reproducible across runs.
%
%   Zero affine distortion: because the reference geometry exactly
%   matches the requested cell, the affine map reduces to identity.
%
% Verified-computation surfaces left untouched:
%   make_mesh_by_gmsh_ref (rigorous gmsh + I_intval mesh construction)
%   affine_map_triangle_mesh (rigorous interval affine transform)
%
% Args:
%   a, b        : triangle parameters (can be intval).
%   h           : mesh size (can be intval).
%   block_size  : kept for backward compatibility; no longer used.
%   verbose     : print cache-miss diagnostic.
%   tag         : optional namespace for the cache key.

    if nargin < 6, tag = ''; end
    if nargin < 5 || isempty(verbose), verbose = false; end
    if nargin < 4, block_size = []; end %#ok<NASGU>  % unused; kept for API stability

    persistent cache
    if isempty(cache)
        cache = containers.Map('KeyType','char','ValueType','any');
    end

    hmid = I_mid(h);
    amid = I_mid(a);
    bmid = I_mid(b);

    % Cache key uses full-precision mid values: a hit <=> identical ref geometry.
    key = sprintf('%s|%.17g|%.17g|%.17g', tag, hmid, amid, bmid);

    if isKey(cache, key)
        st = cache(key);
    else
        st = struct('ref_mesh', [], 'a0', amid, 'b0', bmid);
        st.ref_mesh = make_mesh_by_gmsh_ref(I_intval(amid), I_intval(bmid), h);
        cache(key) = st;
        if verbose
            fprintf('[mesh-cache] new ref mesh key=%s  a0=%.17g  b0=%.17g\n', ...
                key, amid, bmid);
        end
    end

    mesh = affine_map_triangle_mesh(st.ref_mesh, st.a0, st.b0, a, b);
    mesh.domain = I_intval([0 0; 1 0; a b]);
end
