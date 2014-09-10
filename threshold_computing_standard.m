function [upperLimit, innerLimit] = threshold_computing_standard(trainingData)

   
    % Compute mean and std
    m = nanmean(trainingData);
    s = std(trainingData,1);
    
    % Compute limits
    upperLimit = m + 3*s;
    innerLimit = m - 3*s;

end