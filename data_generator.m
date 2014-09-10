% -- Data generator --

% p = Number of sensors

function data = data_generator(p)
    
    % Offset between sensors
    %x_offset = [0 0 0];
    %x_offset = normrnd(0, 10, p, 1);
    %y_offset = [0 5 -5];
    %y_offset = 5*(0:p);
    y_offset = 0*(0:p);

    % Number of instances
    n = 2000;

    % Periods of data
    
    period_active = 1:2000; mu1 = 0; sigma1 = 1;
    period_inactive = {250:450, 1200:1400}; mu2 = -20; sigma2 = 1;
    %period1 = [1 n/4]; mu1 = 20; sigma1 = 1;
    %period2 = [n/4+1 n/2]; mu2 = 0; sigma2 = 1;
    %period3 = [n/2+1 n]; mu3 = 0; sigma3 = 1;

    k = 5;  % shutdown duration
    q = 50; % startup duration
    %shutdown_period = [n/4-k n/4+k];
    %startup_period = [n/2-q n/2+q];

    % Generate data
    data = [];
    for i=1:p
        data(period_active, i) = normrnd(mu1, sigma1, numel(period_active), 1) + y_offset(i);
        for j=1:numel(period_inactive)
            data(period_inactive{j}, i) = normrnd(mu2, sigma2, numel(period_inactive{j}), 1) + y_offset(i);
            
            % Interpolate shutdown
            f = polyfit([period_inactive{j}(1) period_inactive{j}(k)], [mu1 mu2], 1);
            r = period_inactive{j}(1):period_inactive{j}(k);
            data(r, i) = f(1)*r + f(2) + y_offset(i);
            
            %Interpolate startup
            f = polyfit([period_inactive{j}(end-q) period_inactive{j}(end)], [mu2 mu1], 1);
            r = period_inactive{j}(end-q):period_inactive{j}(end);
            data(r, i) = f(1)*r + f(2) + y_offset(i);
        end
        %data(period1(1):period1(2), i) = normrnd(mu1, sigma1, period1(2)-period1(1)+1, 1) + y_offset(i);
        %data(period2(1):period2(2), i) = normrnd(mu2, sigma2, period2(2)-period2(1)+1, 1) + y_offset(i);
        %data(period3(1):period3(2), i) = normrnd(mu3, sigma3, period3(2)-period3(1)+1, 1) + y_offset(i);

        % Fuzzyness of shutdown period
%         for j = shutdown_period(1):shutdown_period(2)
%             if (rand() < 0.5)
%                 data(j,i) = normrnd(mu1, sigma1) + y_offset(i);
%             else
%                 data(j,i) = normrnd(mu2, sigma2) + y_offset(i);
%             end
%         end
        
        % Interpolate
        %f = polyfit([shutdown_period(1) shutdown_period(2)], [mu1 mu2], 1);
        %r = shutdown_period(1):shutdown_period(2);
        %data(r, i) = f(1)*r + f(2) + y_offset(i);

        % Fuzzyness of startup period
%         for j = startup_period(1):startup_period(2)
%             if (rand() < 0.5)
%                 data(j,i) = normrnd(mu2, sigma2) + y_offset(i);
%             else
%                 data(j,i) = normrnd(mu3, sigma3) + y_offset(i);
%             end
%         end
        
        % Interpolate
        %f = polyfit([startup_period(1) startup_period(2)], [mu2 mu3], 1);
        %r = startup_period(1):startup_period(2);
        %data(r, i) = f(1)*r + f(2) + y_offset(i);
    end

    plot(data)
end