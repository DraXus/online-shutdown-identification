% launch_experiments

% Available detection methods:
%  - xie_siegmund_t1
%  - xie_siegmund_t2
%  - mei
%  - tartakovsky_veeravalli
%  - zamba_hawkins
%  - martin-salvador

training_window = 30;
show_plots = 0;
methods = {'XS1', 'XS2', 'MEI', 'TV', 'ZH', 'MS1', 'MS2'};
all_delays = {};
average_delay = zeros(1,numel(methods));
median_delay = zeros(1,numel(methods));

i=1;
for m = methods
    delays = change_detection(m{1}, training_window, show_plots);
    all_delays{i} = delays;
    average_delay(i) = mean(delays);
    median_delay(i) = median(delays);
    i=i+1;
end

num_methods = numel(methods);
f = figure;
group = [];
i=1;
for m = methods
    group = [group repmat({m{1}}, 1, numel(all_delays{i}))];
    i=i+1;
end
boxplot(cell2mat(all_delays), group)