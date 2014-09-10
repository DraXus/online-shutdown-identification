% Asymptotically Optimal Quickest Change Detection in Distributed Sensor Systems
% A. Tartakovsky and V. Veeravalli


function [statistic, change_point, detection_time] = tartakovsky_veeravalli(data, threshold, window_size)
    
    
    % Output values
    statistic = [];
    change_point = NaN;
    detection_time = NaN;

    % Loop control values
    t = 1;
    no_change = 1;
    rows = size(data,1);
    
    while t<rows && no_change
       St = sum(data(1:t, :), 1);
        values = [];
        start = t-window_size;
        if start<1
            start = 1;
        end
        for k=start:t
            Sk = sum(data(1:k, :), 1);
            U = (St - Sk) .* ((t-k)^(-1/2));
            U(isnan(U)) = 0;
            logGLR = (abs(U).^2)/2;
            values = [values sum(logGLR)];
        end
        [maxValue, index] = max(values);
        statistic = [statistic maxValue];

        if maxValue > threshold
           change_point = index;
           detection_time = t;
           no_change = 0;
        end
        
        t = t+1;
    end
    
end
