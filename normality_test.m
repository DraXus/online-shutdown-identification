data_test = data(1:TRAINING_WINDOW,:);
[rows,cols] = size(data_test);

idx = ones(rows,cols);
i=1;
while i<numel(real_changes)
    range = real_changes(i):real_changes(i+1);
    idx(range,:) = 0;
    i = i+2;
end

active = [];
inactive = [];
test = [];
test_inactive = [];
for i=1:cols
    active(:,i) = data_test(idx(:,i)==1,i);
    inactive(:,i) = data_test(idx(:,i)==0,i);
    
    m = nanmean(active(:,i));
    s = nanstd(active(:,i));
    m_inactive = nanmean(inactive(:,i));
    s_inactive = nanstd(inactive(:,i));
    test(i) = kstest((active(:,i)-m)/s);
    test_inactive(i) = kstest((inactive(:,i)-m)/s);
end

x = active(:,3);
[f,x_values] = ecdf(x);
figure;
F = plot(x_values,f);
set(F,'LineWidth',2);
hold on;
G = plot(x_values,normcdf(x_values,0,1),'r-');
set(G,'LineWidth',2);
legend([F G],...
       'Empirical CDF','Standard Normal CDF',...
       'Location','SE');