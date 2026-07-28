function [cells, cover_complete] = upper_build_compact_cells( ...
    y_min, y_max, dx, dy, max_cells)
%UPPER_BUILD_COMPACT_CELLS Rectangular cover of the normalized atlas.
%   For each y-strip we cover
%       1/2 <= x <= sqrt(1-y_min(strip)^2).
%   The final rectangle may extend beyond x^2+y^2=1.  This is intentional:
%   the trial-function inequality remains valid there, so proving the
%   conjecture on this harmless superset closes the curved boundary without
%   special sliver cells.

if ~(0 < y_min && y_min < y_max && y_max < 1 && dx > 0 && dy > 0)
    error('Invalid atlas bounds or cell widths.');
end
if nargin < 5 || isempty(max_cells)
    max_cells = Inf;
end
cover_complete = true;

cells = struct('id', {}, 'x_lo', {}, 'x_hi', {}, 'y_lo', {}, 'y_hi', {});
id = 0;
y_lo = y_min;
while y_lo < y_max
    y_hi = min(y_lo + dy, y_max);
    % Use the outward-rounded upper endpoint itself.  In interval mode this
    % avoids any unproved binary64/ulp assumption at the curved boundary.
    % In double mode I_intval and I_sup reduce to their diagnostic scalar
    % counterparts.
    circle_root = sqrt(I_intval(1)-I_intval(y_lo)^2);
    x_domain_hi = min(1, I_sup(circle_root));
    x_lo = 0.5;
    while x_lo < x_domain_hi
        x_hi = min(x_lo + dx, x_domain_hi);
        id = id + 1;
        cells(id).id = id; %#ok<AGROW>
        cells(id).x_lo = x_lo;
        cells(id).x_hi = x_hi;
        cells(id).y_lo = y_lo;
        cells(id).y_hi = y_hi;
        if id >= max_cells
            cover_complete = false;
            return;
        end
        x_lo = x_hi;
    end
    y_lo = y_hi;
end
end
