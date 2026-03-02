function msh = gmshread(filename)
%GMSHREAD  Fast reader for Gmsh MSH v2.2 ASCII (no Python/meshio).
%
% Supports:
%   - $MeshFormat 2.x, ASCII (filetype=0)
%   - $Nodes (id x y z)
%   - $Elements: type 1 (line), type 2 (triangle)
%
% Returns struct:
%   msh.nodes          (N x 2) double
%   msh.elements       (T x 3) double   triangles
%   msh.edges          (E x 2) double   unique edges of triangles
%   msh.boundary_edges (B x 2) double   boundary edges (from lines or inferred)

    fid = fopen(filename, 'r');
    if fid < 0
        error('gmshread:CannotOpen', 'Cannot open file: %s', filename);
    end
    c = onCleanup(@() fclose(fid)); %#ok<NASGU>

    gotNodes = false;
    gotElems = false;

    node_ids = uint32([]);
    nodes_xy = [];

    tri = zeros(0,3,'uint32');
    bnd = zeros(0,2,'uint32');

    mapMode = 0;      % 0: id==idx, 1: dense array, 2: containers.Map
    id2idx  = uint32([]);
    id2map  = [];

    while true
        line = fgetl(fid);
        if ~ischar(line), break; end
        s = strtrim(line);
        if isempty(s), continue; end

        if strcmp(s,'$MeshFormat')
            fmt = strtrim(fgetl(fid));
            vals = sscanf(fmt,'%f %d %d');
            if numel(vals) < 3
                error('gmshread:BadMeshFormat', 'Bad $MeshFormat line: %s', fmt);
            end
            ver = vals(1);
            filetype = vals(2); % 0=ASCII, 1=binary
            if ver >= 4.0
                error('gmshread:MSH4NotSupported', ...
                      'MSH v%.1f not supported. Output MSH2.2 ASCII (gmsh -format msh2 -bin 0).', ver);
            end
            if filetype ~= 0
                error('gmshread:BinaryNotSupported', ...
                      'Binary MSH2 not supported. Output ASCII (gmsh -format msh2 -bin 0).');
            end
            % consume end tag
            endtag = strtrim(fgetl(fid));
            if ~strcmp(endtag,'$EndMeshFormat')
                % tolerate odd files; try to continue
            end

        elseif strcmp(s,'$Nodes')
            nline = strtrim(fgetl(fid));
            N = sscanf(nline,'%d');
            if isempty(N) || N <= 0
                error('gmshread:BadNodesCount', 'Bad $Nodes count: %s', nline);
            end

            node_ids = zeros(N,1,'uint32');
            nodes_xy = zeros(N,2);

            for i=1:N
                l = fgetl(fid);
                v = sscanf(l,'%d %f %f %f');
                if numel(v) < 4
                    error('gmshread:BadNodeLine','Bad node line: %s', l);
                end
                node_ids(i)  = uint32(v(1));
                nodes_xy(i,1)= v(2);
                nodes_xy(i,2)= v(3);
            end

            endtag = strtrim(fgetl(fid));
            if ~strcmp(endtag,'$EndNodes')
                error('gmshread:BadEndNodes','Expected $EndNodes, got: %s', endtag);
            end
            gotNodes = true;

        elseif strcmp(s,'$Elements')
            if ~gotNodes
                error('gmshread:Order','Encountered $Elements before $Nodes.');
            end

            mline = strtrim(fgetl(fid));
            M = sscanf(mline,'%d');
            if isempty(M) || M < 0
                error('gmshread:BadElementsCount','Bad $Elements count: %s', mline);
            end

            % --- build id -> index map ---
            ids = double(node_ids);
            N = numel(ids);
            if all(ids(:) == (1:N)')
                mapMode = 0;
            else
                maxId = max(ids);
                if maxId <= 5*N
                    mapMode = 1;
                    id2idx = zeros(maxId,1,'uint32');
                    id2idx(node_ids) = uint32(1:N);
                else
                    mapMode = 2;
                    id2map = containers.Map('KeyType','uint32','ValueType','uint32');
                    for k=1:N
                        id2map(node_ids(k)) = uint32(k);
                    end
                end
            end

            tri = zeros(M,3,'uint32'); tri_n = 0;
            bnd = zeros(M,2,'uint32'); bnd_n = 0;

            for e=1:M
                l = fgetl(fid);
                v = sscanf(l,'%d')';
                if numel(v) < 4, continue; end

                etype = v(2);
                ntags = v(3);
                nodeStart = 4 + ntags;

                if etype == 2
                    if numel(v) < nodeStart+2, continue; end
                    n1 = uint32(v(nodeStart));
                    n2 = uint32(v(nodeStart+1));
                    n3 = uint32(v(nodeStart+2));

                    tri_n = tri_n + 1;
                    tri(tri_n,:) = map3(n1,n2,n3);

                elseif etype == 1
                    if numel(v) < nodeStart+1, continue; end
                    n1 = uint32(v(nodeStart));
                    n2 = uint32(v(nodeStart+1));

                    bnd_n = bnd_n + 1;
                    bnd(bnd_n,:) = map2(n1,n2);
                end
            end

            tri = tri(1:tri_n,:);
            bnd = bnd(1:bnd_n,:);

            endtag = strtrim(fgetl(fid));
            if ~strcmp(endtag,'$EndElements')
                error('gmshread:BadEndElements','Expected $EndElements, got: %s', endtag);
            end
            gotElems = true;
        end
    end

    if ~gotNodes || ~gotElems
        error('gmshread:Incomplete','Failed to parse nodes/elements from %s', filename);
    end

    msh = struct();
    msh.nodes = nodes_xy;
    msh.elements = double(tri);

    % boundary edges
    if isempty(bnd)
        E = [tri(:,[1 2]); tri(:,[2 3]); tri(:,[3 1])];
        E = sort(E,2);
        [Eu,~,ic] = unique(E,'rows');
        counts = accumarray(ic,1);
        bnd = Eu(counts==1,:);
    else
        bnd = sort(bnd,2);
        bnd = unique(bnd,'rows');
    end
    msh.boundary_edges = double(bnd);

    % all unique triangle edges
    E = [tri(:,[1 2]); tri(:,[2 3]); tri(:,[3 1])];
    E = sort(E,2);
    msh.edges = double(unique(E,'rows'));

    % ----- local mappers (as nested *anonymous* via closure; no function defs) -----
    function out = map2(i,j) %#ok<DEFNU>
        switch mapMode
            case 0
                out = uint32([i j]);
            case 1
                out = uint32([id2idx(i) id2idx(j)]);
            otherwise
                out = uint32([id2map(i) id2map(j)]);
        end
    end
    function out = map3(i,j,k) %#ok<DEFNU>
        switch mapMode
            case 0
                out = uint32([i j k]);
            case 1
                out = uint32([id2idx(i) id2idx(j) id2idx(k)]);
            otherwise
                out = uint32([id2map(i) id2map(j) id2map(k)]);
        end
    end
end