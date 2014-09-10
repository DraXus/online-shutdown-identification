function [upperLimit, innerLimit] = threshold_computing(trainingData)

    % Compute median and MAD
    m = nanmedian(trainingData);
    s = 1.4826*mad(trainingData,1); %this is actually an estimation of the standard deviation!!
    %s = mad(trainingData,1);
    
    % Compute limits
    upperLimit = m + 3*s;
    innerLimit = m - 3*s;

end