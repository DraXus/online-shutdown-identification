% Sequence multi-sensor change-point detection
% Yao Xie and David Siegmund

% Stopping rule T2

% Assumptions: observations are mutually independent and normally distributed with unit variances.
%              if no change, mean = 0
function [statistic, change_point, detection_time] = xie_siegmund_t2(data, threshold, percentage_affected_sensors, window_size)
    

    data = abs(data);
    
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
            arraylogs = log(1 - percentage_affected_sensors + percentage_affected_sensors*exp(logGLR));
            values = [values sum(arraylogs)];
        end
        [maxValue, index] = max(values);
        statistic = [statistic maxValue];

        if  maxValue > threshold
           change_point = index;
           detection_time = t;
           no_change = 0;
        end
        
        t = t+1;
    end
end