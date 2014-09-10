% -- Multi-sensor change-point detection generic method --
%
% Available detection methods:
%  - XS1: xie_siegmund_t1
%  - XS2: xie_siegmund_t2
%  - MEI: mei
%  - TV: tartakovsky_veeravalli
%  - ZH: zamba_hawkins
%  - MS1: martin_salvador
%  - MS2: martin_salvador2

% Parameters
%  - method: string of the detection method to use
%  - training_window: size of the training window
%  - plots: 0 for no creating plots, 1 otherwise

function delays = change_detection(method, plots)

    %% LOAD DATA
    synthetic = 0;
    random_synthetic = 0;

    if synthetic
        %load('C:\Users\Manuel\BU_data\acrylic_acid\paper_shutdownIdentification\processedData\synthetic\2sensors-2changes.mat', 'data')
        load('### SET INPUT FILE ###', 'data')
        real_changes = [250 450 1200 1400];
        OFFSET = 0;
    elseif random_synthetic
        p = 2; % number of sensors
        data = data_generator(p);
        SENSORS = 1:p;
        OFFSET = 0;
    else
        SENSORS = 5:8;
        %load('C:\Users\Manuel\BU_data\acrylic_acid\paper_shutdownIdentification\processedData\p50-50\rawData1Normalized.mat', 'normalizedData')
        load('### SET INPUT FILE ###', 'normalizedData')
        data = normalizedData(159540:371276, SENSORS);
        %real_changes = [16484 24798 59169 69996 123960 133175 177250 185551];
        real_changes = [123960 133175 177250 185551];
        %training_data = normalizedData(:,SENSORS);
        %OFFSET = size(training_data,1);
        %load('C:\Users\Manuel\BU_data\acrylic_acid\paper_shutdownIdentification\processedData\p50-50\rawData2Normalized.mat', 'normalizedData')
        %test_data = normalizedData(:,SENSORS);
        %load('C:\Users\Manuel\BU_data\acrylic_acid\paper_shutdownIdentification\processedData\p50-50\rawData2NormalizedPeriods.mat', 'auxRealPeriods')
        %real_changes = auxRealPeriods(~isnan(auxRealPeriods));
        %real_changes = real_changes - OFFSET; %offset
        %data = vertcat(training_data, test_data);
        %clear trainig_data test_data normalizedData % free memory
        
        % Sub-sampling real-dataset
        SAMPLING = 10;% 1 sample every 10 minutes
        [rows,cols] = size(data);
        rows = floor(rows/SAMPLING); 
        newdata = zeros(rows, numel(SENSORS));
        for i=1:rows
            newdata(i,:) = data(i*SAMPLING,:);
        end
        data = newdata;
        data = interpolate_missing_values(data);
        clear normalizedData newdata;
        % Update of real_changes
        real_changes = floor(real_changes/SAMPLING); 
    end

    %% DATA PREPROCESSING
    % Interpolate missing values
    data = interpolate_missing_values(data);

    %  Normalize data
    %data = normc(data);

    [rows,cols] = size(data);

    %% DETECTION METHOD
    changes = [];
    detected = [];
    statistic = [];
    begin = 1;
    finish = size(data,1);
    global shutdown;
    shutdown = 0;
    global TRAINING_WINDOW;
    %TRAINING_WINDOW = training_window;
    TRAINING_WINDOW = floor(rows/2);
    
    working_data = data(begin:finish,:);
    % Phase I : Parameter estimation
    %[mean_vector, std_vector] = parameter_estimation(working_data(1:TRAINING_WINDOW, :));
    
    % Phase I : Parameter estimation for each sensor
    distributions = cell(cols,1);
    mean_vector_active = [];
    std_vector_active = [];
    mean_vector_inactive = [];
    std_vector_inactive = [];
    for i=1:cols
        distributions{i} = gmdistribution.fit(working_data(1:TRAINING_WINDOW, i), 2);
        mixing = distributions{i}.PComponents;
        % Data from active periods is more frequent, so probability in the
        % mixing is higher
        if mixing(1) < mixing(2)
            active_dist_index = 2;
            inactive_dist_index = 1;
        else
            active_dist_index = 1;
            inactive_dist_index = 2;
        end
        % Mean and std for active distribution
        mean_vector_active = [mean_vector_active distributions{i}.mu(active_dist_index)];
        std_vector_active = [std_vector_active distributions{i}.Sigma(active_dist_index)];
        % Mean and std for inactive distribution
        mean_vector_inactive = [mean_vector_inactive distributions{i}.mu(inactive_dist_index)];
        std_vector_inactive = [std_vector_inactive distributions{i}.Sigma(inactive_dist_index)];
    end
    
    % Only for martin-salvador method:
    [upperLimit, innerLimit] = threshold_computing(working_data(1:TRAINING_WINDOW, :));
    
    begin = begin+TRAINING_WINDOW;
    statistic = zeros(1, TRAINING_WINDOW);
    
    mean_vector = mean_vector_active;
    std_vector = std_vector_active;
    mean_of_change = mean_vector_inactive;

    %while begin+TRAINING_WINDOW < finish
    while begin < finish
        
%         if ~shutdown
%             mean_vector = mean_vector_active;
%             std_vector = std_vector_active;
%             mean_of_change = mean_vector_inactive;
%         else
%             mean_vector = mean_vector_inactive;
%             std_vector = std_vector_inactive;
%             mean_of_change = mean_vector_active;
%         end
        
        %working_data = normc(data(begin:finish,:));
        working_data = data(begin:finish,:);


        % Phase I : Parameter estimation
        %[mean_vector, std_vector] = parameter_estimation(working_data(1:TRAINING_WINDOW, :));
        % Only for martin-salvador method:
        %[upperLimit, innerLimit] = threshold_computing(working_data(1:TRAINING_WINDOW, :));

        % Phase II : Online detection
        %working_data = working_data(TRAINING_WINDOW+1:end,:);

        try
            switch method

                case 'XS1'
                    % Scale data to zero mean and unit variance
                    %scaled_data = scale_data(working_data, mean_vector, std_vector);
                    %scaled_mean_of_change = scale_data(mean_of_change, mean_vector, std_vector);
                    percentage_affected_sensors = 1;
                    threshold = 500;
                    window_size = 50;
                    %[st, change_point, time] = xie_siegmund_t1(scaled_data, threshold, percentage_affected_sensors, scaled_mean_of_change, window_size);
                    [st, change_point, time] = xie_siegmund_t1(working_data, threshold, percentage_affected_sensors, mean_of_change, window_size);

                case 'XS2'
                    % Scale data to zero mean and unit variance
                    scaled_data = scale_data(working_data, mean_vector, std_vector);
                    percentage_affected_sensors = 1;
                    threshold = 500;
                    window_size = 50;
                    [st, change_point, time] = xie_siegmund_t2(scaled_data, threshold, percentage_affected_sensors, window_size);

                case 'MEI'
                    % Scale data to zero mean and unit variance
                    scaled_data = scale_data(working_data, mean_vector, std_vector);
                    threshold = 500;
                    window_size = 50;
                    [st, change_point, time] = mei(scaled_data, threshold, window_size);

                case 'TV'
                    % Scale data to zero mean and unit variance
                    scaled_data = scale_data(working_data, mean_vector, std_vector);
                    threshold = 500;
                    window_size = 50;
                    [st, change_point, time] = tartakovsky_veeravalli(scaled_data, threshold, window_size);

                case 'ZH'
                    % Scale data to zero mean and unit variance
                    scaled_data = scale_data(working_data, mean_vector, std_vector);
                    alpha = 0.005; % significance level
                    [st, change_point, time] = zamba_hawkins(scaled_data, alpha, mean_vector);

                case 'MS1'
                    warning_window = 30;
                    %working_data = filter(ones(1,50)/50,1,working_data); % filter
                    [st, change_point, time] = onlineShutdownDetectionCUSUM(working_data, upperLimit, innerLimit, warning_window);
                
                case 'MS2'
                    threshold = 0.5;
                    [st, change_point, time] = martin_salvador2(working_data, threshold, distributions);
            end

            %statistic = [statistic zeros(1, TRAINING_WINDOW) st]; %it includes empty statistic for training stage
            statistic = [statistic st];
            if isnan(change_point)
                disp('No more changes found')
                break;
            end
            %changes = [changes begin+TRAINING_WINDOW+change_point-1];
            %detected = [detected begin+TRAINING_WINDOW+time-1];
            %begin = begin + TRAINING_WINDOW + time;
            changes = [changes begin+change_point-1];
            detected = [detected begin+time-1];
            begin = begin+time;         
            if ~shutdown
                shutdown = 1;
            else
                shutdown = 0;
            end
        catch err
            if (strcmp(err.identifier,'MATLAB:EndOfFile'))
              begin = finish;
              disp(err)
           % Display any other errors as usual.
           else
              rethrow(err);
            end
        end
    end



    %% PLOTS
    
    if plots
        fig=figure;
        hax=axes;
        hold on
        plot(statistic(TRAINING_WINDOW+1:end))
        offset_detected = detected - TRAINING_WINDOW;
        title('Statistic')
        for i=1:size(changes,2)
            %line([changes(i) changes(i)],get(hax,'YLim'),'Color',[1 0 0])
            %line([offset_detected(i) offset_detected(i)],get(hax,'YLim'),'Color',[1 0 0])
        end
        hold off

        fig=figure;
        hax=axes;
        hold on
        plot(data(TRAINING_WINDOW+1:end,:))
        title('Data')
        for i=1:size(changes,2)
            %line([changes(i) changes(i)],get(hax,'YLim'),'Color',[1 0 0])
            line([offset_detected(i) offset_detected(i)],get(hax,'YLim'),'Color',[1 0 0])
        end
        hold off
    end
    delays = compute_delays(real_changes, detected);
    average_delay = mean(delays);
    median_delay = median(delays);
    sprintf('Number of changes detected %i (%i uniques)', size(detected,2), size(unique(detected),2))
    sprintf('Average delay %.2f', average_delay)
    sprintf('Median delay %.2f', median_delay)
    sprintf('-------------------------')
    estimated_delays = compute_delays(real_changes, changes);
    average_estimated_delay = mean(estimated_delays);
    median_estimated_delay = median(estimated_delays);
    sprintf('Number of changes detected %i (%i uniques)', size(changes,2), size(unique(changes),2))
    sprintf('Average delay %.2f', average_estimated_delay)
    sprintf('Median delay %.2f', median_estimated_delay)
    sprintf('-------- END ----------')
    
    
end