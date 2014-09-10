function threshold = compute_threshold_xs1(data, mean_vector_inactive, window_size, thresholds, percentage_affected_sensors, num_shutdowns)

    [rows,cols] = size(data);
    statistic_xs1 = [];
    mean_of_change = mean_vector_inactive;
    matrix_mean_of_change = repmat(mean_of_change, window_size, 1);
    matrix_mean_of_change_sq = repmat((mean_of_change.^2)/2, window_size, 1);
    for t=1:rows

        values = [];
        start = t-window_size; % window size
        if start<1
            start = 1;
        end
        for k=start:t
            % Vector of p (# sensors) log-likelihoods of observations accumulated by time t>k
            log_likelihood = sum(matrix_mean_of_change(1:t-k,:).*data(k+1:t,:) - matrix_mean_of_change_sq(1:t-k,:), 1);
            log_likelihood(log_likelihood<0) = 0; % positive part
            % Global log-likelihood of all p sensors
            global_log_likelihood = sum(log(1 - percentage_affected_sensors + percentage_affected_sensors*exp(log_likelihood)));
            values = [values global_log_likelihood];
        end

        maxValue = max(values);
        statistic_xs1 = [statistic_xs1 maxValue];
    end
    

    detected = zeros(1,numel(thresholds));
    shutdown = 0;
    for i=thresholds
        for t=1:rows
            if ~shutdown
                condition = statistic_xs1(t) >= i;
            else
                condition = statistic_xs1(t) < i;
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
    strcat('[XS1] Detected ', num2str(detected(threshold)), ' with threshold=', num2str(threshold))
end