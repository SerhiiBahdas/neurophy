function [] = readaloud(input, tDelay)
%READALOUD Vocalizes the numbers from the input 1-D array or vocalizes 
% strings from the input 1-D string array with a specified pause between 
% the vocalizations. * For Mac users only. 
%
%   INPUT =============================================================
%   
%   input (array)
%   Input array of integers or strings. 
%   Example: [1,2,3,4,5,6] or ["start the treadmill", "stop the treadmill"]
%
%   tDelay (numeric)
%   Time delay between the volcalization of numbers, s. 
%   Example: 1
%
%   AUTHOR ============================================================
%
%   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
%
%   ===================================================================

% Loop through each integer in the vector
for iNumber = 1:length(input)

    % Convert integer to text.
    sNumber = num2str(input(iNumber));
    
    % Create the command to call the macOS 'say' command.
    cmd = sprintf('say %s', sNumber);
    
    % Execute the command.
    system(cmd);
    
    % Optional: Add a pause between numbers.
    pause(tDelay);
end

end % function