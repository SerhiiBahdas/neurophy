function statsStruct = computeStats(relax, bPlot)
%COMPUTESTATS Compute statistical measures for signals in the given structure.
%
%   statsStruct = computeStats(relax, bPlot)
%
%   INPUT ===============================================================
%
%   relax (structure)
%   Structure containing subject-specific fields, each with subfields
%   specific to a muscle, which contain nx1 signals.
%
%   bPlot (boolean)
%   Flag to enable or disable plotting of the signals and statistics.
%
%   OUTPUT ==============================================================
%
%   statsStruct (structure)
%   Structure of the same organization as the input, but containing
%   statistical measures instead of signals.
%
%   ======================================================================
%
%   The function processes the input structure 'relax' and computes the 
%   mean, standard deviation, interquartile range, quartiles, confidence 
%   interval, standard error, and Tukey criterion for each nx1 signal 
%   in the structure.
%
%   If bPlot is set to true, the function will also generate plots for 
%   each signal along with the computed statistical measures.
%
%   Example usage:
%   statsStruct = computeStats(relax, true);
%
%   ======================================================================

% Get the list of subjects from the relax structure
subjects = fieldnames(relax);

% Initialize the output structure
statsStruct = struct();

% Loop over each subject
for iSubj = 1:length(subjects)
    subjName = subjects{iSubj};
    
    % Get the list of muscles for the current subject
    muscles = fieldnames(relax.(subjName));
    
    % Loop over each muscle
    for iMuscle = 1:length(muscles)
        muscleName = muscles{iMuscle};
        
        % Extract the signal for the current subject and muscle
        signal = relax.(subjName).(muscleName);
        
        % Compute the statistical measures
        meanVal = mean(signal);
        stdVals = std(signal) * [1, 2, 3];
        iqrVal = iqr(signal);
        quartiles = quantile(signal, [0.25, 0.5, 0.75]);
        ci = meanVal + tinv([0.025, 0.975], length(signal)-1) * (std(signal)/sqrt(length(signal)));
        stderr = std(signal)/sqrt(length(signal));
        
        % Compute Tukey criterion
        k = 1.5;
        tukeyLower = quartiles(1) - k * iqrVal;
        tukeyUpper = quartiles(3) + k * iqrVal;
        
        % Store the computed measures in the output structure
        statsStruct.(subjName).(muscleName) = struct('Mean', meanVal, ...
                                                     'STD', stdVals, ...
                                                     'IQR', iqrVal, ...
                                                     'Quartiles', quartiles, ...
                                                     'CI', ci, ...
                                                     'StdErr', stderr, ...
                                                     'TukeyLower', tukeyLower, ...
                                                     'TukeyUpper', tukeyUpper);
        
        % Plot the signal and statistics if bPlot is true
        if bPlot
            figure('Position', [10 10 1500 800]);
            plot(signal);
            hold on;
            yline(meanVal, '--', 'Mean');
            yline(meanVal + stdVals(1), '--', '1STD');
            yline(meanVal + stdVals(2), '--', '2STD');
            yline(meanVal + stdVals(3), '--', '3STD');
            yline(quartiles(1), '--', 'Q1');
            yline(quartiles(2), '--', 'Q2');
            yline(quartiles(3), '--', 'Q3');
            yline(ci(1), '--', 'CI Lower');
            yline(ci(2), '--', 'CI Upper');
            yline(tukeyLower, '-.', 'Tukey Lower');
            yline(tukeyUpper, '-.', 'Tukey Upper');
            title([subjName, ' - ', muscleName]);
            hold off;
        end
    end
end

end
