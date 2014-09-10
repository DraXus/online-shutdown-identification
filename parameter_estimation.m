function [mean_vector, std_vector] = parameter_estimation(data)
%     global TRAINING_WINDOW;
    
    %rows = size(data,1);
    
%     if TRAINING_WINDOW > rows
%         throw(MException('MATLAB:EndOfFile', 'Not enough instances to estimate the parameters'));
%     end
    
    mean_vector = mean(data(1:end,:));
    std_vector = std(data(1:end,:));
    
%     % We assume that the change occurs in t > TRAINING_WINDOW
%     if TRAINING_WINDOW+1 > rows
%          % if there is no more instances to process
%         throw(MException('MATLAB:EndOfFile', 'No more instances to process out of the training window'));
%     end
    
end