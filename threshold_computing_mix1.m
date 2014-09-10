function [upperLimit, innerLimit] = threshold_computing_mix1(trainingData)

   
    % Compute mean and std
    m = nanmedian(trainingData);
    s = std(trainingData,1);
    
    % Compute limits
    upperLimit = m + 3*s;
    innerLimit = m - 3*s;

end