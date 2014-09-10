function threshold = compute_threshold_sgz(data, window_size, thresholds, num_shutdowns, upperLimit, innerLimit)

    [rows,cols] = size(data);
    statistic = zeros(1,rows);
    CUSUM_statistic = zeros(window_size,cols);
    variableWeights = ones(1,cols);
    for t=1:rows

        % Outlier?
        is_outlier = isnan(data(t,:)) | data(t,:) < innerLimit | data(t,:) > upperLimit;

        CUSUM_statistic(1:end-1,:) = CUSUM_statistic(2:end,:);
        CUSUM_statistic(end,:) = is_outlier(:);

        median_CUSUM_statistic = median(sum(CUSUM_statistic));
        mad_CUSUM_statistic = mad(sum(CUSUM_statistic),1);
        upper_threshold_3sigma = median_CUSUM_statistic+3*mad_CUSUM_statistic;
        inner_threshold_3sigma = median_CUSUM_statistic-3*mad_CUSUM_statistic;

        disagree = sum(CUSUM_statistic) > upper_threshold_3sigma...
                    | sum(CUSUM_statistic) < inner_threshold_3sigma;

        variableWeights(disagree) = 0;
        variableWeights(~disagree) = 1;

        value = max(sum(CUSUM_statistic(:,variableWeights>0)));
        statistic(t) = value;
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
    
    [~, threshold] = min(abs(detected-num_shutdowns));
    strcat('[SGZ] Detected ', num2str(detected(threshold)), ' with threshold=', num2str(threshold))
end