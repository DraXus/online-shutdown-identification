
% z-score normalization
function data = scale_data(data, mean_vector, std_vector)
    rows = size(data,1);
    data = (data - repmat(mean_vector, rows, 1)).*repmat(1./std_vector, rows, 1);
end

% zero mean
% function data = scale_data(data, mean_vector, std_vector)
%     rows = size(data,1);
%     data = (data - repmat(mean_vector, rows, 1));
% end