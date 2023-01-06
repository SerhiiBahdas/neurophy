function data = extractHeatmapData(varargin)
%EXTRACTHEATMAPDATA Extracts data from the heatmap(s) from the last plotted
% figure. To extract the heatmap data from a specific figure, pass a figure 
% handle to the function. 
%
%   data = extractHeatmapData()
%   data = extractHeatmapData('FigureObj', fig)
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
%   data (cell array)
%   Structure containing X-lables, Y-labels, and color data for all 
%   heatmaps in the figure.
%
%   EXAMPLE ===============================================================
%
%   (Option 1) Extract the data from the last plotted figure:
%
%   data = extractHeatmapData();  
%
%   (Option 2) Pass your figure handle and extract the data: 
%   
%   data = extractHeatmapData('FigureObj', fig); 
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

% Find number of plots in the figure. 
nPlot = length(prop); 

% Loop through plots. 
for iPlot = 1:nPlot 

    % Extract X-lables. 
    data.("Plot" + string(iPlot)).XData = {prop(iPlot).XData}; 
    
    % Extract Y-lables.
    data.("Plot" + string(iPlot)).YData = {prop(iPlot).YData}; 

    % Extarct color data for the heatmap cells.
    data.("Plot" + string(iPlot)).ColorData = {prop(iPlot).ColorData}; 
    
end % iPlot

end % function