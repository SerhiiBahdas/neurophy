function data = extractBoxchartData(varargin)
%EXTRACTBOXCHARTDATA Extracts data from the boxchart(s) from the last plotted
% figure. To extract the boxchart data from a specific figure, pass a figure 
% handle to the function. 
%
%   data = extractBoxchartData()
%   data = extractBoxchartData('FigureObj', fig)
%
%   INPUT =================================================================
%
%   [OPTIONAL] 
% 
%   fig (object)
%   Figure handle. See https://www.mathworks.com/help/matlab/ref/gcf.html 
%   Example: fig = gcf
%
%   OUTPUT ================================================================
%
%   data (structure)
%   Structure containing X-lables and Y-data visualized on all boxchart(s) 
%    in the figure. 
%
%   EXAMPLE ===============================================================
%
%   (Option 1) Extract the data from the last plotted figure:
%
%   data = extractBoxchartData();  
%
%   (Option 2) Pass your figure handle and extract the data: 
%   
%   data = extractBoxchartData('FigureObj', fig); 
%
%   AUTHOR ================================================================
%
%   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
%
%   =======================================================================


%% FETCH ARGUMENTS. Fetch optional arguments. If argument is not specified, 
%   assign it a default value. Check the type of the inputs. 

% By default, fetch the handge of the last-plotted figure. 
fig_default = gcf; 

% Fetch optional arguments. 
p = inputParser; 

% Add optional input argument to a structure and check its type.
addParameter(p, 'FigureObj', fig_default, @isobject);

% Assign the parameters to a structure
parse(p,varargin{:}); p = p.Results;

%% FETCH DATA. Extract graphics object properties

% Query graphics object properties. 
prop = get(p.FigureObj, 'Children'); 

% Extract only objects related to the data (lables are disregarded). 
prop = findobj(prop, 'Type', 'Axes'); 

% Find number of plots in the figure. 
nPlot = length(prop); 

% Loop through plots. 
for iPlot = 1:nPlot 

    % Extract X-lables. 
    data.("Plot" + string(iPlot)).XData = {prop(nPlot+1-iPlot).Children.XData}; 
    
    % Extract Y-data.
    data.("Plot" + string(iPlot)).YData = {prop(nPlot+1-iPlot).Children.YData}; 
    
end % iPlot

end % function