% Interpolate missing values
% Source: http://stackoverflow.com/a/3630003

function new_data = interpolate_missing_values(data)
    
    new_data = [];
    [rows, columns] = size(data);
    for c=1:columns
        time = 1:rows;
        mask = ~isnan(data(:, c));
        new_data(:,c) = data(:,c);
        new_data(~mask,c) = interp1(time(mask), data(mask, c), time(~mask));
    end
end