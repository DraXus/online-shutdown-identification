%GLOBAL_WINDOW_SIZE = 55;

for GLOBAL_WINDOW_SIZE=[55 60 65 70 75 80 85 90 95 100]
%% REAL DATASET
%SENSORS = [4:11 13 15 18];
SENSORS = [4:11 13 15 18];


%load('C:\Users\Manuel\BU_data\acrylic_acid\paper_shutdownIdentification\processedData\p50-50\rawData2Normalized.mat', 'normalizedData')
%data = normalizedData(:, SENSORS);
%clear normalizedData;

%load('C:\Users\Manuel\BU_data\acrylic_acid\paper_shutdownIdentification\rawData\rawDataAllNormed.mat', 'rawDataAllN')
load('### SET INPUT FILE ###', 'rawDataAllN')
data = rawDataAllN(1:1096274, SENSORS);
clear rawDataAllN;

% Adding random noise
%data(:,numel(SENSORS)+1) = rand([1 size(data,1)]);
%data(:,numel(SENSORS)+2) = rand([1 size(data,1)]);

[rows,cols] = size(data);

% Sub-sampling real-dataset
SAMPLING = 10;% 1 sample every 10 minutes
rows = floor(rows/SAMPLING); 
newdata = zeros(rows, cols);

for i=1:rows
    newdata(i,:) = data(i*SAMPLING,:);
end
data = newdata;
clear newdata;

% Interpolate missing values (NaN)
data = interpolate_missing_values(data);
 
TRAINING_WINDOW = floor(rows/2);
%load('C:\Users\Manuel\BU_data\acrylic_acid\paper_shutdownIdentification\processedData\p50-50\rawData2NormalizedPeriods_no-offset.mat', 'realPeriods')
%load('C:\Users\Manuel\BU_data\acrylic_acid\paper_shutdownIdentification\processedData\p50-50\rawData2NormalizedPeriods.mat', 'realPeriods')
load('### SET INPUT FILE ###', 'realPeriods')
real_changes = sort(floor(realPeriods(~isnan(realPeriods))/SAMPLING));
real_changes = real_changes(1:end-1); % removing last change (not-real)
real_changes = real_changes-TRAINING_WINDOW;
clear realPeriods;


%% SYNTHETIC DATASET
%load('C:\Users\Manuel\BU_data\acrylic_acid\paper_shutdownIdentification\processedData\synthetic\2sensors-2changes.mat', 'data')
%[rows,cols] = size(data);
% Simulate failure
%data(800:860,1) = NaN;

%% Phase I : Parameter estimation for each sensor

distributions = cell(cols,1);
mean_vector_active = [];
std_vector_active = [];
mean_vector_inactive = [];
std_vector_inactive = [];
for i=1:cols
    distributions{i} = gmdistribution.fit(data(1:TRAINING_WINDOW, i), 2);
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

%percentage_affected_sensors = 1-min(sum(corr(data(1:TRAINING_WINDOW,:))<0.9))/cols;

percentage_affected_sensors = 1;

%upperLimit = mean_vector_active + 3*std_vector_active;
%innerLimit = mean_vector_active - 3*std_vector_active;

[upperLimit, innerLimit] = threshold_computing(data(1:TRAINING_WINDOW, :)); % for MS1

max_threshold = floor(sum((GLOBAL_WINDOW_SIZE/2)*mean_vector_inactive.^2));
thresholds = 1:max_threshold;
num_shutdowns = 24;
threshold_xs1 = compute_threshold_xs1(data(1:TRAINING_WINDOW, :), mean_vector_inactive, GLOBAL_WINDOW_SIZE, thresholds, percentage_affected_sensors, num_shutdowns);
threshold_xs2 = compute_threshold_xs2(data(1:TRAINING_WINDOW, :), GLOBAL_WINDOW_SIZE, thresholds, percentage_affected_sensors, num_shutdowns);
threshold_mei = compute_threshold_mei(data(1:TRAINING_WINDOW, :), GLOBAL_WINDOW_SIZE, thresholds, num_shutdowns);
threshold_tv = compute_threshold_tv(data(1:TRAINING_WINDOW, :), GLOBAL_WINDOW_SIZE, thresholds, num_shutdowns);

thresholds = 1:GLOBAL_WINDOW_SIZE;
threshold_sgz = compute_threshold_sgz(data(1:TRAINING_WINDOW, :), GLOBAL_WINDOW_SIZE, thresholds, num_shutdowns, upperLimit, innerLimit);

% [upperLimit_hampel, innerLimit_hampel] = threshold_computing_hampel(data(1:TRAINING_WINDOW, :));
% [upperLimit_standard, innerLimit_standard] = threshold_computing_standard(data(1:TRAINING_WINDOW, :));
% [upperLimit_mix1, innerLimit_mix1] = threshold_computing_mix1(data(1:TRAINING_WINDOW, :));
% [upperLimit_mix2, innerLimit_mix2] = threshold_computing_mix2(data(1:TRAINING_WINDOW, :));


%weights = sum(corr(data(1:TRAINING_WINDOW,:)))/cols; % for MS2

data = data(TRAINING_WINDOW+1:end,:);
[rows,cols] = size(data);

% Simulate sensor failure
%data(11000:12000, 1) = mean_vector_inactive(1);

delays = {};
methods = {'XS1', 'XS2', 'MEI', 'TV', 'ZH', 'SGZ', 'MS2'};

%% XS1
disp('XS1')
%profile on
statistic_xs1 = [];
change_points_xs1 = [];
detected_xs1 = [];
mean_of_change = mean_vector_inactive;
%percentage_affected_sensors = 1;
shutdown = 0;
threshold = threshold_xs1;
window_size = GLOBAL_WINDOW_SIZE;
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
        %log_likelihood = sum(repmat(mean_of_change, t-k, 1).*data(k+1:t,:) - repmat((mean_of_change.^2)/2, t-k, 1), 1);
        log_likelihood = sum(matrix_mean_of_change(1:t-k,:).*data(k+1:t,:) - matrix_mean_of_change_sq(1:t-k,:), 1);
        log_likelihood(log_likelihood<0) = 0; % positive part
        % Global log-likelihood of all p sensors
        global_log_likelihood = sum(log(1 - percentage_affected_sensors + percentage_affected_sensors*exp(log_likelihood)));
        values = [values global_log_likelihood];
    end

    [maxValue, index] = max(values);
    statistic_xs1 = [statistic_xs1 maxValue];
    
    if ~shutdown
        condition = maxValue >= threshold;
    else
        condition = maxValue < threshold;
    end
    
    if condition
        %change_points_xs1 = [change_points_xs1 index+start];
        detected_xs1 = [detected_xs1 t];
        shutdown = ~shutdown;
    end
    
end
delays{1} = compute_delays(real_changes, detected_xs1);
%profile viewer
%profile off
%load gong.mat;
%soundsc(y(1:10000))


%% XS2
disp('XS2')
%profile on
statistic_xs2 = [];
change_points_xs2 = [];
detected_xs2 = [];
%percentage_affected_sensors = 1;
shutdown = 0;
threshold = threshold_xs2;
window_size = GLOBAL_WINDOW_SIZE;
for t=1:rows
    %St = sum(data(1:t, :), 1);
    values = [];
    start = t-window_size;
    if start<1
        start = 1;
    end
    
    for k=start:t
        %Sk = sum(data(1:k, :), 1);
        %U = (St - Sk) .* ((t-k)^(-1/2));
        %U(isnan(U)) = 0;
        if t>k
            U = sum(data(k+1:t, :), 1)/(t-k);
            %logGLR = (U(U>=0).^2)/2;
            logGLR = (abs(U).^2)/2;
            arraylogs = log(1 - percentage_affected_sensors + percentage_affected_sensors*exp(logGLR));
            values = [values sum(arraylogs)];
        end
    end
    [maxValue, index] = max(values);
    statistic_xs2 = [statistic_xs2 maxValue];
    
    if ~shutdown
        condition = maxValue >= threshold;
    else
        condition = maxValue < threshold;
    end
    
    if condition
       % change_points_xs2 = [change_points_xs2 index+start];
        detected_xs2 = [detected_xs2 t];
        shutdown = ~shutdown;
    end
end
delays{2} = compute_delays(real_changes, detected_xs2);
%profile viewer
%profile off
%load gong.mat;
%soundsc(y(1:10000))

%% MEI
disp('MEI')
%profile on
statistic_mei = [];
change_points_mei = [];
detected_mei = [];
shutdown = 0;
threshold = threshold_mei;
window_size = GLOBAL_WINDOW_SIZE;
for t=1:rows
    %St = sum(data(1:t, :), 1);
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
    [maxValue, index] = max(values); % maximum by each column
    st = sum(maxValue);
    statistic_mei = [statistic_mei st];
    
    if ~shutdown
        condition = st >= threshold;
    else
        condition = st < threshold;
    end
    
    if condition
        %change_points_mei = [change_points_mei index];
        detected_mei = [detected_mei t];
        shutdown = ~shutdown;
    end
end
delays{3} = compute_delays(real_changes, detected_mei);
% profile viewer
% profile off
% load gong.mat;
% soundsc(y(1:10000))

%% TV
disp('TV')
% profile on
statistic_tv = [];
change_points_tv = [];
detected_tv = [];
shutdown = 0;
threshold = threshold_tv;
window_size = GLOBAL_WINDOW_SIZE;
for t=1:rows
    %St = sum(data(1:t, :), 1);
    values = [];
    start = t-window_size;
    if start<1
        start = 1;
    end
    for k=start:t
        if t>k
            %Sk = sum(data(1:k, :), 1);
            %U = (St - Sk) .* ((t-k)^(-1/2));
            %U(isnan(U)) = 0;
            U = sum(data(k+1:t, :), 1)/(t-k);
            logGLR = (abs(U).^2)/2;
            values = [values sum(logGLR)];
        end
    end
    [maxValue, index] = max(values);
    statistic_tv = [statistic_tv maxValue];
    
    if ~shutdown
        condition = maxValue >= threshold;
    else
        condition = maxValue < threshold;
    end
    
    if condition
        %change_points_tv = [change_points_tv index];
        detected_tv = [detected_tv t];
        shutdown = ~shutdown;
    end
end
delays{4} = compute_delays(real_changes, detected_tv);
% profile viewer
% profile off
% load gong.mat;
% soundsc(y(1:10000))

%% ZH
% disp('ZH')
% profile on
% alpha = 0.005; % significance level
% statistic_zh = [];
% %change_points_zh = [];
% detected_zh = [];
% shutdown = 0;
% threshold = 10;
% window_size = 25;
% for t=1:rows
%     values = [];
%     start = t-window_size;
%     if start<1
%         start = 1;
%     end
%     for k=start:t
%         if t>k+1
%             m1k = mean(data(start:k,:));
%             mkt = mean(data(k+1:t,:));
%             yk = (m1k - mkt) * sqrt(k*(t-k)/t);
%             
%             w1k = 0;
%             for i=start:k
%                 w1k = w1k + (data(i,:)-m1k)*(data(i,:)-m1k)';
%             end
% 
%             wkt = 0;
%             for i=k+1:t
%                 wkt = wkt + (data(i,:)-mkt)*(data(i,:)-mkt)';
%             end
% 
%             wk = (w1k+wkt)/(t-2);
%             t2st = (yk/wk)*yk';
% 
%             values = [values t2st];
%         end
%     end
%     [maxValue, index] = max(values);
%     statistic_zh = [statistic_zh maxValue];
%     
%     if ~shutdown
%         condition = maxValue >= threshold;
%     else
%         condition = maxValue < threshold;
%     end
%     
%     if condition
%         detected_zh = [detected_zh t];
%         shutdown = ~shutdown;
%     end
% end
% delays{5} = compute_delays(real_changes, detected_zh);
% profile viewer
% profile off
delays{5} = [];
%load gong.mat;
%soundsc(y(1:10000))

%% SGZ
disp('SGZ')
% profile on
statistic_ms = [];
%maxCUSUM = 30;
%maxConsecutiveAnomaly = maxCUSUM/2;
window_size = GLOBAL_WINDOW_SIZE;
CUSUM_statistic = zeros(window_size,cols);
consecutiveInlier = zeros(1,cols);
variableWeights = ones(1,cols);
historicalWeights = zeros(rows,cols);
change_points_ms1 = [];
detected_ms1 = [];
shutdown = 0;
threshold = threshold_sgz;
for t=1:rows
    % Outlier?
    is_outlier = isnan(data(t,:)) | data(t,:) < innerLimit | data(t,:) > upperLimit;
    % Inlier?
    is_inlier = ~is_outlier;
    
    %consecutiveInlier = consecutiveInlier.*is_inlier + is_inlier;
    %idx = consecutiveInlier < maxConsecutiveAnomaly;
    %CUSUM_statistic(~idx) = 0;
    %CUSUM_statistic(idx) = CUSUM_statistic(idx) + is_outlier(idx);
    
    %CUSUM_statistic = CUSUM_statistic.*is_outlier + is_outlier;
    %for i=1:cols
    CUSUM_statistic(1:end-1,:) = CUSUM_statistic(2:end,:);
    CUSUM_statistic(end,:) = is_outlier(:);
    %end
    
    % Dynamic adjusting of weights
%     median_CUSUM_statistic = median(CUSUM_statistic);
%     mad_CUSUM_statistic = mad(CUSUM_statistic,1)./0.6745;
%     upper_threshold_3sigma = median_CUSUM_statistic+3*mad_CUSUM_statistic;
%     inner_threshold_3sigma = median_CUSUM_statistic-3*mad_CUSUM_statistic;
     
    median_CUSUM_statistic = median(sum(CUSUM_statistic));
    mad_CUSUM_statistic = mad(sum(CUSUM_statistic),1);
    upper_threshold_3sigma = median_CUSUM_statistic+3*mad_CUSUM_statistic;
    inner_threshold_3sigma = median_CUSUM_statistic-3*mad_CUSUM_statistic;
    
    % If disagree with majority -> disable variable
%     disagree = CUSUM_statistic > upper_threshold_3sigma...
%                 | CUSUM_statistic < inner_threshold_3sigma;
            
    disagree = sum(CUSUM_statistic) > upper_threshold_3sigma...
                | sum(CUSUM_statistic) < inner_threshold_3sigma;

    variableWeights(disagree) = 0;
    variableWeights(~disagree) = 1;
    
    historicalWeights(t,:) = variableWeights;
    
    %value = max(CUSUM_statistic(variableWeights>0));
    value = max(sum(CUSUM_statistic(:,variableWeights>0)));
    
    if ~shutdown
        condition = value >= threshold;
    else
        condition = value < threshold;
    end
    
    statistic_ms = [statistic_ms value];

    if condition
        change_points_ms1 = [change_points_ms1 t];
        detected_ms1 = [detected_ms1 t];
        shutdown = ~shutdown;
    end
end
delays{6} = compute_delays(real_changes, detected_ms1);
% profile viewer
% profile off
% load gong.mat;
% soundsc(y(1:10000))

%% MS2
% disp('MS2')
% statistic_ms2 = {};
% for i=1:cols
%     statistic_ms2{i} = [];
% end
% statistic_ms2_active = [];
% change_points_ms2 = [];
% detected_ms2 = [];
% shutdown = 0;
% threshold = 0.8;
% for t=1:rows
%     p_active = zeros(1,cols);
%     p_inactive = zeros(1,cols);
%     for i=1:cols
%         if ~isnan(data(t,i))
%             p_both = posterior(distributions{i}, data(t,i));
% 
%             mixing = distributions{i}.PComponents;
%             % Data from active periods is more frequent, so probability in the
%             % mixing is higher
%             if mixing(1) < mixing(2)
%                 active_dist_index = 2;
%             else
%                 active_dist_index = 1;
%             end
% 
%             p_active(i) = p_both(active_dist_index);
%             statistic_ms2{i} = [statistic_ms2{i} p_active(i)];
%         end
%     end
%     st = sum(p_active)/cols;
%     statistic_ms2_active = [statistic_ms2_active st];
%     
%     if ~shutdown
%         condition = st < threshold; 
%     else
%         condition = st >= threshold;
%     end
% 
%     if condition
%         change_points_ms2 = [change_points_ms2 t];
%         detected_ms2 = [detected_ms2 t];
%         shutdown = ~shutdown;
%     end
% end
% delays{7} = compute_delays(real_changes, detected_ms2);
delays{7} = [];
% load gong.mat;
% soundsc(y(1:10000))

%% Boxplots of delays
% num_methods = numel(methods);
% f = figure;
% group = [];
% i=1;
% for m = methods
%     group = [group repmat({m{1}}, 1, numel(delays{i}))];
%     i=i+1;
% end
% boxplot(cell2mat(delays), group)
% 
% par_delays = []; par_group = [];
% impar_delays = []; impar_group = [];
% all_delays = cell2mat(delays);
% for i=1:numel(all_delays)
%     if mod(i,2)==0
%         par_delays = [par_delays all_delays(i)];
%         par_group = [par_group group(i)];
%     else
%         impar_delays = [impar_delays all_delays(i)];
%         impar_group = [impar_group group(i)];
%     end
% end
% fig = figure;
% subplot(2,1,1); boxplot(impar_delays, impar_group); title('Shutdowns delays')
% subplot(2,1,2); boxplot(par_delays, par_group); title('Startups delays')

%% Save data
%save(strcat('C:\Users\Manuel\Desktop\experiments\results-',num2str(GLOBAL_WINDOW_SIZE),'.mat'),...
%    'delays', 'statistic_*', 'GLOBAL_WINDOW_SIZE', 'threshold_*', 'detected_*');
save('### SET OUTPUT FILE ###',...
    'delays', 'statistic_*', 'GLOBAL_WINDOW_SIZE', 'threshold_*', 'detected_*');

%% Sorting boxplots for the paper
% TV(4) MEI(3) XS1(1) XS2(2) SGZ(6)
group = [];
all_delays = [];
idx_methods = [4 3 1 2 6];

for i=idx_methods
    group = [group repmat({methods{i}}, 1, numel(delays{i}))];
    all_delays = [all_delays delays{i}];
end

par_delays = []; par_group = [];
impar_delays = []; impar_group = [];
for i=1:numel(all_delays)
    if mod(i,2)==0
        par_delays = [par_delays all_delays(i)];
        par_group = [par_group group(i)];
    else
        impar_delays = [impar_delays all_delays(i)];
        impar_group = [impar_group group(i)];
    end
end

figure1 = figure('Color',[1 1 1]);
% Create axes
axes1 = axes('Parent',figure1,'XTickLabel','','XTick',zeros(1,0),...
    'Position',[0.0703971119133574 0.551799200507035 0.913357400722022 0.407721635912239],...
    'FontSize',12,...
    'FontName','Times New Roman');
hggroup1 = hggroup('Parent',axes1);
boxplot(impar_delays, impar_group, 'Parent', hggroup1);
title('Shutdowns delays','FontSize',12,'FontName','Times New Roman');

axes2 = axes('Parent',figure1,'XTickLabel','','XTick',zeros(1,0),...
    'Position',[0.0703971119133574 0.0414110429447854 0.913357400722022 0.422021686403196],...
    'FontSize',12,...
    'FontName','Times New Roman');
hggroup2 = hggroup('Parent',axes2);
boxplot(par_delays, par_group, 'Parent', hggroup2);
title('Startups delays','FontSize',12,'FontName','Times New Roman');

%print(figure1, strcat('C:\Users\Manuel\Desktop\experiments\plot-',num2str(GLOBAL_WINDOW_SIZE),'.fig'));
print(figure1, '### SET OUTPUT FILE ###');
end