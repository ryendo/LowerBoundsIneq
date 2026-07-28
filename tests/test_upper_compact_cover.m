function test_upper_compact_cover()
%TEST_UPPER_COMPACT_COVER Check rigorous containment of the circular edge.
project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(project_root,'src','upper_conjecture'));
addpath(fullfile(project_root,'src','interval'));

if ~strcmp(getenv('RUN_INTLAB_SMOKE'),'1')
    fprintf(['upper compact-cover interval test skipped; set ', ...
        'RUN_INTLAB_SMOKE=1 and INTLAB_ROOT to enable it.\n']);
    return;
end

upper_prepare_runtime(project_root,'interval');
[cells,complete] = upper_build_compact_cells(0.84,0.851,0.002,0.003,Inf);
assert(complete);
assert(~isempty(cells));

y_lowers = unique([cells.y_lo]);
for k = 1:numel(y_lowers)
    y_lo = y_lowers(k);
    strip = cells([cells.y_lo] == y_lo);
    [~,order] = sort([strip.x_lo]);
    strip = strip(order);

    assert(strip(1).x_lo == 0.5);
    for j = 1:(numel(strip)-1)
        assert(strip(j).x_hi == strip(j+1).x_lo);
    end

    exact_circle = sqrt(I_intval(1)-I_intval(y_lo)^2);
    assert(strip(end).x_hi >= I_sup(exact_circle));
end

fprintf('test_upper_compact_cover: PASS (%d cells, %d strips)\n', ...
    numel(cells),numel(y_lowers));
end
