function threshold = compute_threshold_mei(data, window_size, thresholds, num_shutdowns)

    [rows,cols] = size(data);
    statistic = zeros(1,rows);
    for t=1:rows

        values = [];
        start = t-window_size;
        if start<1
            start = 1;
        end
        for k=start:t
            if t>k
                U = sum(data(k+1:t, :), 1)/(t-k);
                %Sk = sum(data(1:k, :), 1);
                %U = (St - Sk) .* ((t-k)^(-1/2));
                %U(isnan(U)) = 0;
                logGLR = (abs(U).^2)/2;
                values = [values 
                          logGLR];
            end
        end
        
        if numel(values)>0
            [maxValue, index] = max(values); % maximum by each column
            st = sum(maxValue);
            statistic(t) = st;
        end
    end
    

    detected = zeros(1,numel(thresholds));
    shutdown = 0;
    for i=thresholds
        for t=1:rows
            if ~shutdown
                condition = statistic(t) >= i;
            else
                condition = statistic(t) < i;
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
    strcat('[MEI] Detected ', num2str(detected(threshold)), ' with threshold=', num2str(threshold))
end