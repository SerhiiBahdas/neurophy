function sFileList = getContent(sDir,varargin)
%GETCONTENT Lists all folders and/or files in the specified directory. The
%   extension of the files can be specified. 
%
%   sFileList = getContent(sDir)
%   sFileList = getContent(sDir, 'Extension', sFileExt)
%
%   INPUT ===========================================================
%
%   sDir (string)
%   The path to the folder with contents. 
%   Example: "C:\Users\sb0220\Desktop"
%
%   sExtension (char, string)
%   Specify the extension of the files you want to list in the chosen 
%       directory. By default, function lists all folders and files.    
%   Example: getContent(..., 'Extension', '*.mat'); 
%
%   OUTPUT ==========================================================
%
%   sFileList (string)
%   List of the folders or/and files in the chosen directory. 
%
%   EXAMPLE =========================================================
%
%   Get all folders and files in the current directory:
%
%   sDir = "C:\Users\sb0220\Desktop"; 
%   sFileList = getContent(sDir); 
%
%   Get files with specific file extension: 
%   
%   sFileExt = '*.mat'; 
%   sDir = "C:\Users\sb0220\Desktop"; 
%   sFileList = getContent(sDir, 'Extension', sFileExt);
%
%   AUTHOR ==========================================================
%
%   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
%
%   =================================================================

%% DEFAULT VALUES. By default, the file extension is set to 'all contents'
%   (including folders and files).

sDefaultExt = '*'; 

%% FETCH ARGUMENTS. Fetch optional arguments. If argument is not specified, 
%   assign it a default value. Check the type of the inputs. 

% Fetch optional arguments
p = inputParser; 

% Function handle testing if the input argument is either char or string
% array
isword = @(sWord) ischar(sWord) || isstring(sWord); 

% Add required argument to a structure
addRequired(p,'sDir',isword);

% Add optional input argument to a structure and check its type
addParameter(p, 'Extension', sDefaultExt, isword);

% Assign the parameters to a structure
parse(p,sDir,varargin{:});

%% LIST CONTENTS. Assign folder contents to the output string array.

% Fetch the folders or/and files with specified extension
contents = dir(fullfile(p.Results.sDir, p.Results.Extension));

% Remove the rows containing '.' and '..' in the name column
contents = contents(~ismember({contents(:).name}, {'.','..', '.DS_Store'}));

% Create a string array containing the names of the folders or/and
%   files in in the specified directory 
sFileList = convertCharsToStrings(extractfield(contents,'name'))';

end % function

