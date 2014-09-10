
% assuming that both real_changes and detected_changes are sorted
function delays = compute_delays(real_changes, detected_changes)

    delays = [];
    
    i = 1; % index of real_changes
    j = 1; % index of detected_changes
    num_real_changes = numel(real_changes);
    num_detected_changes = numel(detected_changes);
    
    while j <= num_detected_changes

        % case 1
        if detected_changes(j) <= real_changes(i)
            negative_delay = detected_changes(j) - real_changes(i);
            delays = [delays negative_delay];
            j = j+1;
        
        % case 2
        elseif detected_changes(j) > real_changes(i) &&...
               i+1 <= num_real_changes &&...
               detected_changes(j) >= real_changes(i+1)
            i = i+1;
       
        % case 3
        elseif detected_changes(j) > real_changes(i) &&...
               i+1 <= num_real_changes &&...
               detected_changes(j) < real_changes(i+1)
            positive_delay = detected_changes(j) - real_changes(i);
            negative_delay = detected_changes(j) - real_changes(i+1);
            if abs(negative_delay) < positive_delay
                delays = [delays negative_delay];
            else
                delays = [delays positive_delay];
            end
            j = j+1;
            
        % case 4
        else
            positive_delay = detected_changes(j) - real_changes(i);
            delays = [delays positive_delay];
            j = j+1;
        end

    end
end