function threshold = compute_threshold_xs2(data, window_size, thresholds, percentage_affected_sensors, num_shutdowns)

    [rows,cols] = size(data);
    statistic_xs2 = zeros(1,rows);
    for t=1:rows

        values = [];
        start = t-window_size;
        if start<1
            start = 1;
        end

        for k=start:t
            if t>k
                U = sum(data(k+1:t, :), 1)/(t-k);
                logGLR = (abs(U).^2)/2;
                arraylogs = log(1 - percentage_affected_sensors + percentage_affected_sensors*exp(logGLR));
                values = [values sum(arraylogs)];
            end
        end
        if numel(values)>0
            [maxValue, index] = max(values);
            statistic_xs2(t) = maxValue;
        end
    end
    

    detected = zeros(1,numel(thresholds));
    shutdown = 0;
    for i=thresholds
        for t=1:rows
            if ~shutdown
                condition = statistic_xs2(t) >= i;
            else
                condition = statistic_xs2(t) < i;
            end

            if condition
                detected(i) = detected(i)+1;
                shutdown = ~shutdown;
            end
        end
        if detected(i) == num_shutdowns
            break;
        end
    end
    
    [minValue, threshold] = min(abs(detected-num_shutdowns));
    strcat('[XS2] Detected ', num2str(detected(threshold)), ' with threshold=', num2str(threshold))
end