function [lo, hi] = upper_interval_endpoints(value, mode)
%UPPER_INTERVAL_ENDPOINTS Convert a scalar enclosure to binary64 endpoints.
if strcmpi(mode, 'interval')
    lo = inf(value);
    hi = sup(value);
else
    lo = value;
    hi = value;
end
lo = double(lo);
hi = double(hi);
end
