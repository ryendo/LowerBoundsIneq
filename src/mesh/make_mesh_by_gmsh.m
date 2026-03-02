function mesh = make_mesh_by_gmsh(a, b, h, block_size, verbose, tag)
% Every block_size calls, generate a reference mesh on (0,0)-(1,0)-(a0,b0),
% then reuse it via affine mapping for subsequent calls.

    if nargin < 6, tag = ''; end
    if nargin < 5 || isempty(verbose), verbose = false; end
    if nargin < 4 || isempty(block_size), block_size = 50; end

    persistent cache
    if isempty(cache)
        cache = containers.Map('KeyType','char','ValueType','any');
    end

    hmid = I_mid(h);
    key = sprintf('%s|%.16g', tag, hmid);

    if ~isKey(cache, key)
        st = struct('counter', 0, 'ref_mesh', [], 'a0', NaN, 'b0', NaN);
        cache(key) = st;
    end
    st = cache(key);

    st.counter = st.counter + 1;

    need_new = isempty(st.ref_mesh) || mod(st.counter-1, block_size) == 0;

    if ~need_new && ~isnan(st.b0)
        beta_now = I_mid(b) / st.b0;
        if beta_now > 1.2 || beta_now < 1/1.2
            need_new = true;
        end
    end

    if need_new
        st.a0 = I_mid(a);
        st.b0 = I_mid(b);
        
        st.ref_mesh = make_mesh_by_gmsh_ref(I_intval(st.a0), I_intval(st.b0), h);

        if verbose
            fprintf('[mesh-cache] new ref mesh key=%s  cnt=%d  a0=%.17g  b0=%.17g\n', ...
                key, st.counter, st.a0, st.b0);
        end
    end
    
    cache(key) = st;
    mesh = affine_map_triangle_mesh(st.ref_mesh, st.a0, st.b0, a, b);
    mesh.domain = I_intval([0 0; 1 0; a b]);    
end