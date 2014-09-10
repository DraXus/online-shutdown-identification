% Scale data between 0 and 1
function data = minmax_scale(data)
    rows = size(data,1);
    data = (data - repmat(min(data),rows,1)).*(1./(repmat(max(data), rows, 1) - repmat(min(data), rows, 1)));
end