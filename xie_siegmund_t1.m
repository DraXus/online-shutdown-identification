% Sequence multi-sensor change-point detection
% Yao Xie and David Siegmund

% Stopping rule T1

% Assumptions: observations are mutually independent and normally distributed with unit variances.
%              if no change, mean = 0
function [statistic, change_point, detection_time] = xie_siegmund_t1(data, threshold, percentage_affected_sensors, mean_of_change, window_size)
    
    % This method only works with positives changes of the mean
    % Therefore, we get the absolute value of data
    data = abs(data);
    mean_of_change = abs(mean_of_change);
    
    % Output values
    statistic = [];
    change_point = NaN;
    detection_time = NaN;
    
    % Loop control values
    t = 1;
    no_change = 1;
    rows = size(data,1);
    
    global shutdown;
    
    while t<rows && no_change
        
        values = [];
        start = t-window_size; % window size
        if start<1
            start = 1;
        end
        for k=start:t
            % Vector of p (# sensors) log-likelihoods of observations accumulated by time t>k
            log_likelihood = sum(repmat(mean_of_change, t-k, 1).*data(k+1:t,:) - repmat((mean_of_change.^2)/2, t-k, 1), 1);
            log_likelihood(log_likelihood<0) = 0; % positive part
            % Global log-likelihood of all p sensors
            arraylogs = sum(log(1 - percentage_affected_sensors + percentage_affected_sensors*exp(log_likelihood)));
            values = [values arraylogs];
        end
        
        [maxValue, index] = max(values);
        statistic = [statistic maxValue];

        if ~shutdown
            condition = maxValue > threshold;
        else
            condition = maxValue <= threshold;
        end
        
        if condition
           change_point = index;
           detection_time = t;
           no_change = 0;
        end
        
        t = t+1;
        
    end
