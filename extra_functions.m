%data = normalizedData(1:200000,5);
data = data_generator(1);
window = 50;
rows = size(data,1);
newdata = zeros(1,rows);

% Sliding window mean
for i=1:rows
    if(i<=window)
        newdata(i) = data(i);
    else
        newdata(i) = mean(data(i-window:i));
    end
end

% Cumulative mean
for i=1:rows
    newdata(i) = mean(data(1:i));
end

% Log-likelihood Gaussian Mixture Model
obj = gmdistribution.fit(data,2);
m0 = obj.mu(1);
m1 = obj.mu(2);
s0 = obj.Sigma(1);
s1 = obj.Sigma(2);
S = zeros(1,rows);
for i=1:rows
    p0 = pdf('Normal', data(i), m0, s0);
    p1 = pdf('Normal', data(i), m1, s1);
    S(i) = log(p0/p1);
    if(i<=window)
        newdata(i) = S(i);
    else
        newdata(i) = sum(S(i-window:i));
    end
end

% Log-likelihood (windowed)
S = zeros(1,rows);
m0 = median(data);
s0 = mad(data);
for i=1:rows
    
    if(i<=window)
        newdata(i) = 0; % equiprobable
    else
        %m0 = mean(data(i-window:i));
        %s0 = std(data(i-window:i));
        p0 = pdf('Normal', data(i), m0, s0);
        S(i) = log(p0/(1-p0));
        newdata(i) = sum(S(i-window:i));
    end
end