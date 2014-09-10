% Manuel Martin Salvador

function [statistic, change_point, detection_time] = onlineShutdownDetectionCUSUM(data, upperLimit, innerLimit, warning_window)
    
    
    [rows, cols] = size(data);
    maxCUSUM = warning_window;
    maxConsecutiveAnomaly = maxCUSUM/2;
    DYNAMIC_WEIGHTS = 1;

   
    CUSUM_statistic = zeros(1,cols);
    consecutiveInlier = zeros(1,cols);
    %consecutiveOutlier = zeros(1,cols);
    variableWeights = ones(1,cols);
    

    % Output values
    statistic = [];
    change_point = NaN;
    detection_time = NaN;


    % Loop control values
    t = 1;
    no_change = 1;
    rows = size(data,1);
    global shutdown;
    
    %% DETECTION LOOP

    while t<rows && no_change
            instance = data(t,:);
            %instance = faultSimulation(instance);
            
            % Outlier?
            is_outlier = isnan(instance) | instance < innerLimit | instance > upperLimit;
            % Inlier?
            is_inlier = ~is_outlier;
            
            if shutdown
                is_outlier = ~is_outlier;
                is_inlier = ~is_inlier;
            end
            
            % if outlier -> increment counter
            consecutiveInlier = consecutiveInlier.*is_inlier + is_inlier;
            idx = consecutiveInlier < maxConsecutiveAnomaly;
            CUSUM_statistic(~idx) = 0;
            CUSUM_statistic(idx) = CUSUM_statistic(idx) + is_outlier(idx);
            
            if DYNAMIC_WEIGHTS
                % Dynamic adjusting of weights
                median_CUSUM_statistic = median(CUSUM_statistic);
                mad_CUSUM_statistic = mad(CUSUM_statistic,1)./0.6745;
                upper_threshold_3sigma = median_CUSUM_statistic+3*mad_CUSUM_statistic;
                inner_threshold_3sigma = median_CUSUM_statistic-3*mad_CUSUM_statistic;
                %upper_threshold_2sigma = median_CUSUM_statistic+2*mad_CUSUM_statistic;
                %inner_threshold_2sigma = median_CUSUM_statistic-2*mad_CUSUM_statistic;

                % If disagree with majority -> disable variable
                disagree = CUSUM_statistic > upper_threshold_3sigma...
                            | CUSUM_statistic < inner_threshold_3sigma;
                variableWeights(disagree) = 0;
                variableWeights(~disagree) = 1;

                % If agree with majority -> enable variable
%                 agree = ~disagree & (CUSUM_statistic < upper_threshold_2sigma...
%                             & CUSUM_statistic > inner_threshold_2sigma);
%                 variableWeights(agree) = 1;

            end
            
            if sum(variableWeights)<1
                %WTF!
                throw('No relevant sensors!! :(')
            end
            
            
            if ~shutdown
                % For detecting shutdowns, we prefer a fast approach for
                % detecting them as soon as possible
                value = max(CUSUM_statistic(variableWeights>0));
                condition = value >= maxCUSUM;
            else
                % For detecting startups, we prefer to wait until the whole
                % process is working.
                value = min(CUSUM_statistic(variableWeights>0));
                condition = value >= maxCUSUM*2; % startups take longer
            end
            
            statistic = [statistic value];

            if condition
                change_point = t-warning_window;
                detection_time = t;
                no_change = 0;
            end

            % Timestamp increase
            t = t + 1;
    end

end