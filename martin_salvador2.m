function [statistic, change_point, detection_time] = martin_salvador2(data, threshold, distributions)

    % Output values
    statistic = [];
    change_point = NaN;
    detection_time = NaN;


    % Loop control values
    t = 1;
    no_change = 1;
    [rows,cols] = size(data);
    
    active_dist_index = [];
    for i=1:cols
        mixing = distributions{i}.PComponents;
        % Data from active periods is more frequent, so probability in the
        % mixing is higher
        if mixing(1) < mixing(2)
            active_dist_index(i) = 2;
        else
            active_dist_index(i) = 1;
        end
    end
    
    global shutdown;
    
    while t<rows && no_change
        
        p_active = zeros(1,cols);
        
        for i=1:cols
            if ~isnan(data(t,i))
                p_both = posterior(distributions{i}, data(t,i));
                p_active(i) = p_both(active_dist_index(i));
            else
                p_active(i) = 0;
            end
        end
        value = sum(p_active)/cols;
        if shutdown
            value = 1-value;
        end
        statistic = [statistic value];
        
        if value < threshold
           change_point = t;
           detection_time = t;
           no_change = 0;
        end
        t = t+1;
        
   end

end

