%% Statistics of MS2
fig = figure;
ax = [];
ax(1) = subplot(cols+1,1,1); plot(data); title('DATA', 'FontWeight', 'bold', 'FontSize', 11)
colors = get(gca,'ColorOrder');
number_of_colors = size(colors,1);
for i=1:cols
    ax(i+1) = subplot(cols+1,1,i+1); 
    if i>number_of_colors
        c = mod(i,number_of_colors);
    else
        c = i;
    end
    plot(statistic_ms2{i}, 'Color', colors(c,:), 'LineWidth', 1.5)
    title(sprintf('Sensor %i', i), 'FontWeight', 'bold', 'FontSize', 11)
end
linkaxes(ax,'x');
axis tight


%% Statistics of methods
fig = figure;
ax = [];
num_methods = 6;
ax(1) = subplot(num_methods,1,1); plot(data); title('DATA', 'FontWeight', 'bold', 'FontSize', 11)
ax(2) = subplot(num_methods,1,2); plot(statistic_mei, 'LineWidth', 1.5); title('MEI', 'FontWeight', 'bold', 'FontSize', 11)
ax(3) = subplot(num_methods,1,3); plot(statistic_tv, 'LineWidth', 1.5); title('TV', 'FontWeight', 'bold', 'FontSize', 11)
ax(4) = subplot(num_methods,1,4); plot(statistic_xs1, 'LineWidth', 1.5); title('XS1', 'FontWeight', 'bold', 'FontSize', 11)
ax(5) = subplot(num_methods,1,5); plot(statistic_xs2, 'LineWidth', 1.5); title('XS2', 'FontWeight', 'bold', 'FontSize', 11)

%ax(6) = subplot(8,1,6); plot(statistic_zh, 'LineWidth', 1.5); title('ZH', 'FontWeight', 'bold', 'FontSize', 11)
ax(6) = subplot(num_methods,1,6); plot(statistic_ms, 'LineWidth', 1.5); title('SGZ', 'FontWeight', 'bold', 'FontSize', 11)
%ax(8) = subplot(8,1,8); plot(statistic_ms2_active, 'LineWidth', 1.5); title('MS2', 'FontWeight', 'bold', 'FontSize', 11)
linkaxes(ax,'x');
axis tight

%% Historical weights
fig = figure;
ax = [];
plot_cols = ceil(sqrt(cols+1));

ax(1) = subplot(plot_cols,plot_cols,1:plot_cols); plot(data); title('DATA', 'FontWeight', 'bold', 'FontSize', 11)
colors = get(gca,'ColorOrder');
number_of_colors = size(colors,1);
for i=1:cols
    ax(i+1) = subplot(plot_cols,plot_cols,i+plot_cols);
    if i>number_of_colors
        c = mod(i,number_of_colors);
    else
        c = i;
    end
    area(historicalWeights(:,i), 'FaceColor', colors(c,:))
    title(sprintf('Sensor %i', i), 'FontWeight', 'bold', 'FontSize', 11)
end
linkaxes(ax,'x');
axis tight


%% Correlation over time
k=500;
correlations = [];
for t=1:rows
    if t<=k
        c = ones(1,cols);
    else
        c = sum(corr(data(t-k:t,:)))/cols;
    end
    correlations = [correlations 
                    c];
end
fig = figure;
ax = [];
for i=1:cols
    ax(i) = subplot(cols, 1, i);
    plot(correlations(:,i))
end
linkaxes(ax,'x');
axis tight

%% Thresholds
sensor = 1;
range = 1:TRAINING_WINDOW;
fig = figure;
ax = [];
rows = range(end);

ax(1) = subplot(3,2,1);
hold on
title('median +- 3*1.4826*MAD')
plot(data(range,sensor), 'Color', [0.5 0.5 0.5])
line([1 rows], [upperLimit(sensor) upperLimit(sensor)])
line([1 rows], [innerLimit(sensor) innerLimit(sensor)])
hold off

ax(2) = subplot(3,2,2);
hold on
title('median +- 3*MAD')
plot(data(range,sensor), 'Color', [0.5 0.5 0.5])
line([1 rows], [upperLimit_hampel(sensor) upperLimit_hampel(sensor)])
line([1 rows], [innerLimit_hampel(sensor) innerLimit_hampel(sensor)])
hold off

ax(3) = subplot(3,2,3);
hold on
title('mean +- 3*std')
plot(data(range,sensor), 'Color', [0.5 0.5 0.5])
line([1 rows], [upperLimit_standard(sensor) upperLimit_standard(sensor)])
line([1 rows], [innerLimit_standard(sensor) innerLimit_standard(sensor)])
hold off

ax(4) = subplot(3,2,4);
hold on
title('median +- 3*std')
plot(data(range,sensor), 'Color', [0.5 0.5 0.5])
line([1 rows], [upperLimit_mix1(sensor) upperLimit_mix1(sensor)])
line([1 rows], [innerLimit_mix1(sensor) innerLimit_mix1(sensor)])
hold off

ax(5) = subplot(3,2,5);
hold on
title('mean +- 3*MAD')
plot(data(range,sensor), 'Color', [0.5 0.5 0.5])
line([1 rows], [upperLimit_mix2(sensor) upperLimit_mix2(sensor)])
line([1 rows], [innerLimit_mix2(sensor) innerLimit_mix2(sensor)])
hold off

ax(6) = subplot(3,2,6);
hold on
title('mean +- 3*std (from gmdistribution)')
upperLimit_gm = mean_vector_active + 3*std_vector_active;
innerLimit_gm = mean_vector_active - 3*std_vector_active;
plot(data(range,sensor), 'Color', [0.5 0.5 0.5])
line([1 rows], [upperLimit_gm(sensor) upperLimit_gm(sensor)])
line([1 rows], [innerLimit_gm(sensor) innerLimit_gm(sensor)])
hold off

linkaxes(ax,'x');
axis tight


%% Two thresholds
sensor = 4;
range = 1:TRAINING_WINDOW;
range_test = TRAINING_WINDOW+1:size(data,1);
fig = figure;
ax = [];
rows = range(end);

ax(1) = subplot(2,2,1);
hold on
title('TRAINING ---- median +- 3*1.4826*MAD')
plot(data(range,sensor), 'Color', [0.5 0.5 0.5])
line([1 rows], [upperLimit(sensor) upperLimit(sensor)])
line([1 rows], [innerLimit(sensor) innerLimit(sensor)])
hold off
ax(2) = subplot(2,2,2);
hold on
title('TRAINING ---- mean +- 3*std (from gmdistribution)')
upperLimit_gm = mean_vector_active(sensor) + 3*std_vector_active(sensor);
innerLimit_gm = mean_vector_active(sensor) - 3*std_vector_active(sensor);
plot(data(range,sensor), 'Color', [0.5 0.5 0.5])
line([1 rows], [upperLimit_gm upperLimit_gm])
line([1 rows], [innerLimit_gm innerLimit_gm])
hold off
ax(1) = subplot(2,2,3);
hold on
title('TEST ---- median +- 3*1.4826*MAD')
plot(data(range_test,sensor), 'Color', [0.5 0.5 0.5])
line([1 rows], [upperLimit(sensor) upperLimit(sensor)])
line([1 rows], [innerLimit(sensor) innerLimit(sensor)])
hold off
ax(2) = subplot(2,2,4);
hold on
title('TEST ---- mean +- 3*std (from gmdistribution)')
upperLimit_gm = mean_vector_active(sensor) + 3*std_vector_active(sensor);
innerLimit_gm = mean_vector_active(sensor) - 3*std_vector_active(sensor);
plot(data(range_test,sensor), 'Color', [0.5 0.5 0.5])
line([1 rows], [upperLimit_gm upperLimit_gm])
line([1 rows], [innerLimit_gm innerLimit_gm])
hold off

linkaxes(ax,'x');
axis tight


%% REAL CHANGES
fig=figure;
hax=axes;
hold on
plot(data(TRAINING_WINDOW+1:end,:))
title('Data')
for i=1:numel(real_changes)
    %line([changes(i) changes(i)],get(hax,'YLim'),'Color',[1 0 0])
    line([real_changes(i) real_changes(i)],get(hax,'YLim'),'Color',[1 0 0])
end
hold off