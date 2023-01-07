function data = extractCurveData(varargin)
%EXTRACTCURVEDATA Extracts data from the 2-D plot(s) in the figure.
%
%   data = extractCurveData()
%   data = extractCurveData('FigureObj', fig, 'bMultPlot', 0)
%
%   INPUT =================================================================
%
%   [OPTIONAL] 
% 
%   FigureObj (object)
%   Figure handle. See https://www.mathworks.com/help/matlab/ref/gcf.html
%   By default function works with the last opened figure. 
%   Example: fig = gcf
%
%   bMultPlot (boolean)
%   If the figure contains multiple plots, set to 1. Otherwise set to 0.
%   The default value is 0. 
%   Example: 0 
%
%   OUTPUT ================================================================
%
%   data (structure)
%   Structure containing data from the figure. 
%
%   EXAMPLE ===============================================================
%
%   (Option 1) Extract data from the last plotted figure:
%
%   data = extractCurveData();  
%
%   (Option 2) Pass your figure handle to a specific figure: 
%   
%   data = extractCurveData('FigureObj', fig); 
%
%   (Option 3) Pass 1 if your figure contains subplots:
%
%   data = extractCurveData('FigureObj', fig, 'bMultPlot', 1)
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

% Create a function that checks if the value is boolean.
isbool = @(x) x ==1 || x==0; 

% By default, assume figure contains only one plot. 
bMultPlot_default = 0; 

% Fetch optional arguments. 
p = inputParser; 

% Add optional input argument to a structure and check its type.
addParameter(p, 'FigureObj', fig_default, @isobject);

% Add optional input argument to a structure and check its type.
addParameter(p, 'bMultPlot', bMultPlot_default, isbool);

% Assign the parameters to a structure
parse(p,varargin{:}); p = p.Results;

%% FETCH DATA. Extract graphics object properties

% Query graphics object properties. 
prop = get(p.FigureObj, 'Children'); 

% If there are multiple plots.
if p.bMultPlot == 1

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

% If there is only one plot.
else 

    % Extract only objects related to the data (lables are disregarded). 
    prop = findobj(prop, 'Type', 'Line');

    % Find number of plots in the figure. 
    nObj = length(prop); 
    
    % Loop through objects. 
    for iObj = 1:nObj 
    
        % Extract X-lables. 
        data.Plot1.XData(iObj) = {prop(nObj+1-iObj).XData}; 
        
        % Extract Y-data.
        data.Plot1.YData(iObj) = {prop(nObj+1-iObj).YData}; 
        
    end % iPlot

end % bMultPlot

end % function