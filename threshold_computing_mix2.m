function [upperLimit, innerLimit] = threshold_computing_mix2(trainingData)

   
    % Compute mean and std
    m = nanmean(trainingData);
    s = mad(trainingData,1);
    
    % Compute limits
    upperLimit = m + 3*s;
    innerLimit = m - 3*s;

end