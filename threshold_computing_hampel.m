function [upperLimit, innerLimit] = threshold_computing_hampel(trainingData)

    % Compute median and MAD
    m = nanmedian(trainingData);
    %s = mad(trainingData,1)./0.6745; %this is actually the standard deviation!!
    s = mad(trainingData,1);
    
    % Compute limits
    upperLimit = m + 3*s;
    innerLimit = m - 3*s;

end