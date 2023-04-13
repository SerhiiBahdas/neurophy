function dependencies = findDependencies(filename)
% FINDDEPENDENCIES Finds the dependencies of a MATLAB script on other files.
%
%   dependencies = findDependencies(filename)
%
%   INPUT =================================================================
%   
%   filename (char array)
%   The name of the MATLAB script file (including extension) to analyze.
%   Example: [1,2,3,4,5,6]
%
%   OUTPUT =================================================================
%
%   dependencies (cell array)
%   Full file paths that the input script depends on.
%
%   AUTHOR ================================================================
%
%   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
%
%   =======================================================================


% Validate input
if ~ischar(filename) || isempty(filename)
    error('Input must be a non-empty character array representing the file name.');
end

% Check if the file exists
if ~exist(filename, 'file')
    error('The specified file does not exist.');
end

% Get the dependencies using the requiredFilesAndProducts function.
try
    [dependencies, ~] = matlab.codetools.requiredFilesAndProducts(filename);
catch ME
    error('Error while analyzing dependencies: %s', ME.message);
end

end % function
