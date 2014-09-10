all_windows = [20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100];
all_methods = [4 3 1 2 6]; % TV(4) MEI(3) XS1(1) XS2(2) SGZ(6)
shutdown_median_delays = cell(numel(all_windows), numel(all_methods));
startup_median_delays = cell(numel(all_windows), numel(all_methods));
false_alarms = cell(numel(all_windows), numel(all_methods));
thresholds = cell(numel(all_windows), numel(all_methods));
num_changes = 22;

idx_window = 1;
for GLOBAL_WINDOW_SIZE=all_windows
    
    %load(strcat('C:\Users\Manuel\Desktop\experiments\results-',num2str(GLOBAL_WINDOW_SIZE),'.mat'))
    load('### SET OUTPUT FILE ###')
    idx_method = 1;
    all_detections = [numel(detected_tv) numel(detected_mei) numel(detected_xs1) numel(detected_xs2) numel(detected_ms1)];
    all_thresholds = [threshold_tv threshold_mei threshold_xs1 threshold_xs2 threshold_sgz];
    for i=all_methods
       par_delays = [];
       impar_delays = [];
       for j=1:numel(delays{i})
           if mod(j,2)==0
               par_delays = [par_delays delays{i}(j)];
           else
               impar_delays = [impar_delays delays{i}(j)];
           end
       end
       aux1 = median(par_delays);
       aux2 = median(impar_delays);
       strcat('[', num2str(i), '] Median startups: ', num2str(aux1),...
           ' Median shutdowns: ', num2str(aux2))
       shutdown_median_delays{idx_window,idx_method} = aux2;
       startup_median_delays{idx_window,idx_method} = aux1;
       
       false_alarms{idx_window,idx_method} = all_detections(idx_method)-num_changes;
       thresholds{idx_window,idx_method} = all_thresholds(idx_method);
       idx_method = idx_method+1;
    end
    
    idx_window = idx_window+1;
    
end

