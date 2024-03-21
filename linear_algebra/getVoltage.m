%GETVOLTAGE Simulates the voltage distribution in a 2D space.
%
%   [V] = getVoltage('gridLen', value, 'gridHeight', value, 'numPointX',...
%        value, 'numPointY', value, 'I', value, 'R', value, 'wireStart',...
%        [x, y], 'wireEnd', [x, y], 'tol', value, 'bPlot', value)
%
%   INPUT =================================================================
%
%   gridLen (numeric)
%   The physical length of the simulation grid in meters.
%   Example: 1
%
%   gridHeight (numeric)
%   The physical height of the simulation grid in meters.
%   Example: 1
%
%   numPointX (numeric)
%   The number of grid points along the length.
%   Example: 50
%
%   numPointY (numeric)
%   The number of grid points along the height.
%   Example: 50
%
%   I (numeric)
%   The current through the wire in Amperes.
%   Example: 1
%
%   R (numeric)
%   The resistance of the wire in Ohms.
%   Example: 1
%
%   wireStart (array)
%   The start point of the wire in meters [x, y].
%   Example: [0.25, 0.4]
%
%   wireEnd (array)
%   The end point of the wire in meters [x, y].
%   Example: [0.25, 0.6]
%
%   tol (numeric)
%   The tolerance for convergence based on changes in voltage.
%   Example: 1e-5
%
%   bPlot (logical)
%   Flag to plot the results. Set to true to plot, false to not.
%   Example: true
%
%   OUTPUT ============================================================
%
%   V (matrix)
%   The 2D array representing the voltage distribution around the wire.
%
%   PLOTTING ==========================================================
%
%   If 'bPlot' is set to true, generates a heat map and a contour plot 
%   showing the voltage distribution around the wire, including a visual 
%   representation of the wire's location.
%
%   EXAMPLE ============================================================
%
%   % Simulate and plot the voltage distribution in a 1m x 1m area with a 
%   % wire running from [0.25, 0.4] to [0.25, 0.6] meters, using a grid of
%   % 50x50 points, with 1A current, 1 Ohm resistance, and plotting the 
%   % results.
%
%   V = getVoltage('gridLen', 1, 'gridHeight', 1, 'numPointX', 50,...
%                  'numPointY', 50, 'I', 1, 'R', 1,...
%                  'wireStart', [0.25, 0.4], 'wireEnd', [0.25, 0.6],...
%                  'tol', 1e-5, 'bPlot', true);
%
%   AUTHOR ============================================================
%
%   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
%
%   ===================================================================


function V = getVoltage(varargin)
    

%% PARAMETER INITIALIZATION. Use input parser to fetch inputs. 

    % Initialize input parser
    p = inputParser;
    
    % Default values
    defLen    = 1;           % meters
    defHeight = 1;           % meters
    defGridX  = 50;         % points
    defGridY  = 50;         % points
    defI      = 1;           % Amperes
    defR      = 1;           % Ohms
    defStart  = [0.25, 0.4]; % meters
    defEnd    = [0.25, 0.6]; % meters
    defTol    = 1e-3;
    defPlot   = true;
    
    % Adding parameters
    addParameter(p, 'gridLen',    defLen,    @isnumeric);
    addParameter(p, 'gridHeight', defHeight, @isnumeric);
    addParameter(p, 'numPointX',  defGridX,  @isnumeric);
    addParameter(p, 'numPointY',  defGridY,  @isnumeric);
    addParameter(p, 'I',          defI,      @isnumeric);
    addParameter(p, 'R',          defR,      @isnumeric);
    addParameter(p, 'wireStart',  defStart,  @isnumeric);
    addParameter(p, 'wireEnd',    defEnd,    @isnumeric);
    addParameter(p, 'tol',        defTol,    @isnumeric);
    addParameter(p, 'bPlot',      defPlot,   @islogical);
    
    % Parsing inputs
    parse(p, varargin{:});
    
    % Extracting results
    gridLen    = p.Results.gridLen;
    gridHeight = p.Results.gridHeight;
    numPointX  = p.Results.numPointX;
    numPointY  = p.Results.numPointY;
    I          = p.Results.I;
    R          = p.Results.R;
    wireStart  = p.Results.wireStart;
    wireEnd    = p.Results.wireEnd; 
    tol        = p.Results.tol;
    bPlot      = p.Results.bPlot;

    %% 
    
    % Convert wireStart and wireEnd from meters to grid points
    nScaleX    = gridLen/numPointX; 
    nScaleY    = gridHeight/numPointY; 

    startPt = [wireStart(1)/nScaleX, wireStart(2)/nScaleY];
    endPt   = [wireEnd(1)/nScaleX, wireEnd(2)/nScaleY];

    % Initialize voltage matrix
    V = zeros(numPointX, numPointY);
    
    % Set initial voltage condition
    Vini = I * R;
    V(startPt(1):endPt(1), startPt(2):endPt(2)) = Vini;
    
    % Iterative update
    deltaV = 1;
    while deltaV > tol
        Vnew = V;
        for i = 1:numPointX
            for j = 1:numPointY
                if ~(i >= startPt(1) && i <= endPt(1) && j >= startPt(2) && j <= endPt(2)) % Check if outside wire area
                    % Determine the indices for the neighboring points, mirroring at boundaries
                    iLeft  = i - 1; if iLeft < 1, iLeft = i + 1; end
                    iRight = i + 1; if iRight > numPointX, iRight = i - 1; end
                    jDown  = j - 1; if jDown < 1, jDown = j + 1; end
                    jUp    = j + 1; if jUp > numPointY, jUp = j - 1; end
                    
                    % Average the voltages of the four neighbors, using mirrored values at boundaries
                    Vnew(i, j) = mean([V(iRight, j), V(iLeft, j), V(i, jUp), V(i, jDown)]);
                end
            end
        end
        
        % Calculate the maximum change in voltages across the grid
        deltaV = max(max(abs(V - Vnew)));
        V = Vnew;
    end

    
% Plotting
if bPlot
    figure;
    hold on;

    % Define the grid points for X and Y based on physical dimensions and number of points
    x = linspace(0, gridLen, numPointX);
    y = linspace(0, gridHeight, numPointY);
    
    % Voltage heatmap
    imagesc(x, y, V');
    colormap(flip(sky)); % Using 'hot' for intensity visualization
    colorbar;
    axis xy; % Ensures the y-axis starts from the bottom
    
    % Contour lines for voltage levels, using the x and y vectors for correct dimensionality
    [C, h] = contour(x, y, V', 'y'); % Yellow contours for visibility
    clabel(C, h); % Labeling contour lines with voltage values
    
    % Wire location on plot, directly using physical dimensions
    plot([wireStart(1), wireEnd(1)], [wireStart(2), wireEnd(2)], 'w-', 'LineWidth', 2); % Wire in white for visibility

    title('Voltage Distribution with Contour Lines');
    xlabel('Length (m)');
    ylabel('Height (m)');

    hold off;
end



