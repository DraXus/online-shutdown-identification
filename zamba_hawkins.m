% A multivariate change-point model for statistical process control
% K.D. Zamba and D.M. Hawkins

function [statistic, change_point, detection_time] = zamba_hawkins(data, alpha, mean)

    [rows,cols] = size(data);
    
    % Output values
    statistic = zeros(1,cols);
    change_point = NaN;
    detection_time = NaN;

    % Loop control values
    t = 1;
    no_change = 1;
    
    
    % T2 Hotelling test
    while t<rows && no_change
        if t>cols
            [T2,P] = T2Hot1(data(1:t,:), alpha, mean);
            statistic = [statistic P];
            if P < alpha
                change_point = t;
                detection_time = t;
                no_change = 0;
            end
        end
        t = t+1;
    end
    
end