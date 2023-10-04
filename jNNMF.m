function [nProj, jBasis] = jNNMF(nData, varargin)
%JNNMF Perform j-Derivative Non-Negative Matrix Factorization on input data.
%
%   [nProj, jBasis] = JNNMF(nData, Name, Value) performs j-Non-Negative
%   Matrix Factorization on nData and returns the projections (nProj) and 
%   j-Basis vectors (jBasis). It also provides visualization options and 
%   allows for the inclusion of additional parameters through name-value 
%   pair arguments.
%
%   INPUT =================================================================
%
%   nData (numeric matrix)
%       The input data matrix with non-negative entries that is to be 
%           factored. 
%       Each row is a variable and each column an observation.
%       Example: randn(10,1000)
%
%   Name-Value Pair Arguments
%   Specify optional comma-separated pairs of Name,Value arguments. Name is
%   the argument name and Value is the corresponding value. 
%
%       - 'numComponents'  : Number of components to be computed 
%                           (DEFAULT = 4)
%       - 'bPlanes'        : Logical flag for visualizing 2D projections on 
%                            all combinations of two jVec planes 
%                            (DEFAULT = true)
%       - 'bAnimate'       : Logical flag for visualizing comet plots 
%                            (DEFAULT = true)
%       - 'b3D'            : Logical flag for performing 3D visualization 
%                            is greater than or equal to 3 
%                            (DEFAULT = false)
%
%   OUTPUT ================================================================
%
%   nProj (numeric matrix)
%       The projected data onto the derived j-Basis vectors.
%
%   jBasis (numeric matrix)
%       The derived j-Basis vectors OF j-DNNMF.
%
%   FUNCTION OPERATIONS ===================================================
%
%   (1) Validate and parse input parameters.
%   (2) Ensure the data is non-negative by adding absolute of minimum negative
%       values to respective columns.
%   (3) Execute NNMF and derive j-Basis vectors.
%   (4) Optionally visualize results based on input flags:
%           a. 3D visualizations of projections on first 3 j-Basis vectors.
%           b. 2D comet plot visualizations.
%           c. 2D scatter plots on all combinations of two j-Basis planes.
%
%   EXAMPLE ===============================================================
% 
%    % DATA MATRIX PARAMETERS.
%    numDim     = 40;    % number of dimensions
%    numSamples = 10000; % number of samples
% 
%    % SIGNAL PARAMETERS.
%    nA     = 1;     % semi-major axis
%    nB     = 0.5;   % semi-minor axis
%    nFreq  = 10;    % frequency
%    nPhase = pi/4;  % phase
% 
%    % GENERATE 2-D TRAJECTORY.
%    tTime = linspace(0, 5*pi, numSamples);
%    nTraj  = [nA * cos(nFreq * tTime); nB * sin(nFreq * tTime + nPhase)];
% 
%    % PROJECT SIGNALS TO 10-D.
%    nB = orth(randn(numDim, 2));
%    nX = (nB * nTraj);
%
%   [nProj, jBasis] = jNNMF(nX)
%
%   AUTHOR =================================================================
%
%   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
%
%   =======================================================================


%% DEFAULT VALUES. Assign default values here. 
paramConfig = {
    'numComponents', 4, @isnumeric;
    'bPlanes',       1, @(x) islogical(x) || x==1 || x==0;
    'bAnimate',      1, @(x) islogical(x) || x==1 || x==0;
    'b3D',           0, @(x) islogical(x) || x==1 || x==0;
};

%% FETCH ARGUMENTS. Fetch optional arguments and check their type. If 
% argument is not specified, assign a default value.  

inputParams = inputParser;

% Add parameters to the parser via a loop.
for idxVec = 1:size(paramConfig, 1)
    addParameter(inputParams, paramConfig{idxVec,1}, ...
        paramConfig{idxVec,2}, paramConfig{idxVec,3});
end

% Assign the parsed parameters to variables.
parse(inputParams, varargin{:});
paramResults = inputParams.Results;

% Direct assignment for better readability in code usage.
numComponents   = paramResults.numComponents;     
bPlanes         = paramResults.bPlanes;     
bAnimate        = paramResults.bAnimate;     
b3D             = paramResults.b3D;     

% Request even number of principal components. 
if mod(numComponents,2)
    error('Specify even number of principal components (numPC).')
end % mod

%% PERFORM NNMF and further calculations.

% Set the random seed
rng(10);

% To ensure nData is non-negative, find the minimum value in the data matrix.
nMinValList = min(nData);

% Find the indices of the negative minimum values.
idxNegativeMinVals = nMinValList < 0;

% Add the negative minimum values to corresponding columns of nData.
nData(:, idxNegativeMinVals) = nData(:, idxNegativeMinVals) + ...
    abs(nMinValList(idxNegativeMinVals));

% Performing Non-Negative Matrix Factorization.
[~, H] = nnmf(nData, numComponents);

% Calculate derivative of the coefficient matrix.
dH = diff(H, 1, 2)';

% Adjust dimensions: Remove the last column from coefficient matrix.
H = H(:, 1:end-1)';

% Calculate the least squares solution M.
M = pinv(H' * H) * H' * dH; 

% Calculate and display skew-symmetric part of M.
Mskew = (M - M') / 2;

disp('Skew-symmetric part of M:'); disp(Mskew);

% Optimize skew-symmetric transformation.
Mskew_optimized = optimSkew(Mskew, H, dH);

disp('Optimized skew-symmetric matrix:'); 
disp(Mskew_optimized);

% Perform eigen-decomposition of Mskew.
[V, D] = eig(Mskew_optimized);

% Order eigenvectors by the imaginary part of eigenvalues.
[~, idxOrder] = sort(imag(diag(D)), 'descend');

V = V(:, idxOrder);

% Initialize jFactor vectors.
jBasis = zeros(size(H, 2));

% Calculate real-valued jFactor vectors from complex eigenvectors.
for idxVec = 1:2:numComponents
    jBasis(:, idxVec)   = real(V(:, idxVec)) - imag(V(:, idxVec+1));
    jBasis(:, idxVec+1) = real(V(:, idxVec)) + imag(V(:, idxVec+1));
end

% Project the score matrix onto the jVec space.
nProj = H * jBasis;


%% VISUALIZE. Visualization Based on Flags

% (1) 3D Visualization
if b3D == 1 && numComponents >= 3
    visualize3DProjections(nProj);
end

% (2) Comet plots for all combinations of nProj
if bAnimate == 1
    createCometPlots(nProj);
end

% (3) Visualize projections on all combinations of two jVec planes
if bPlanes == 1
    visualizePlanesProjections(nProj);
end

end % jNNMF

%% ========================================================================
%                  _   _                      
%  _ __ ___   __ _| |_| |__                   
% | '_ ` _ \ / _` | __| '_ \                  
% | | | | | | (_| | |_| | | |                 
% |_| |_| |_|\__,_|\__|_| |_|                 
%   __                  _   _                 
%  / _|_   _ _ __   ___| |_(_) ___  _ __  ___ 
% | |_| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
% |  _| |_| | | | | (__| |_| | (_) | | | \__ \
% |_|  \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
       
function M_opt = optimSkew(M, X, Xdot, nLearningRate, numIter, nEps)
    %OPTIMSKEW Optimize symmetric skew matrix using gradient descent.
    %
    %   INPUT =======================================================
    %
    %   M (numeric array)
    %   Initial guess for a symmetric skew matrix. 
    %   Example: M = rand(2); M = (M - M') / 2; 
    %
    %   X (numeric array)
    %   Data matrix. 
    %   Example: X = rand(100,2); 
    %
    %   Xdot (numeric array)
    %   Differential of the data matrix. 
    %   Xdot = diff(X); X = X(end-1,:); 
    %
    %   [OPTIONAL PARAMETERS]
    %
    %   nLearningRate (numeric)
    %   Learning rate for gradient descent. 
    %   Example: 1e-4
    %
    %   numIter (nuemric)
    %   Maximum number of iterations for gradient descent. 
    %   Example: 1000
    %
    %   nEps (numeric)
    %   Eror tollerance for gradient descent method. 
    %   Example: 1e-3
    %
    %   OUTPUT ======================================================
    %
    %   M_opt (numeric array)
    %   Optimized symmetric skew matrix. 
    %
    %
    %   AUTHOR ======================================================
    %
    %   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
    %
    %   =============================================================

    % Validate dimensions of input matrices X and Xdot.
    if size(X, 1) ~= size(Xdot, 1)
        error('The number of rows in X and Xdot must be equal');
    end

    % Initialize parameters for gradient descent
    if nargin < 4
        % Set default learning rate.
        nLearningRate = 1e-3; 
    end

    if nargin < 5
        % Set default number of iterations.
        numIter = 1000; 
    end

    if nargin < 6
        % Set default tolerance for change in cost function.
        nEps = 1e-6; 
    end

    % Define the cost function and a gradient.
    function [J, G] = cost(M, X, Xdot)
        X_M = X * M;

        % Cost function. 
        J = norm(Xdot - X_M, 'fro');

        % Unconstrained gradient.
        G_raw = computePartialDerivative(Xdot, X, M);

        % Project the gradient onto the space of skew-symmetric matrices.
        G = 0.5 * (G_raw - G_raw'); 

    end % cost
    
    % Intialize a variable storing output of the cost function. 
    nCost_previous = Inf;

    disp('Optimizing the symmetric skew matrix...')

    % Run gradient descent algorithm.
    for iItter = 1:numIter

        % Every 100 iterations. 
        if ~mod(iItter,100)

            % Display every 100-th iteration. 
            disp("iteration..." + string(iItter))

            % Reduce learning rate. 
            nLearningRate = nLearningRate/10; 

        end

        % If the maximum itteration reached. 
        if iItter == numIter
            disp('Maximum number of iterations was reached.')
        end

        % Compute cost and gradient. 
        [nCost_current, G] = cost(M, X, Xdot);

        % Update matrix. 
        M = M - nLearningRate * G;
        
        % Check if change in cost function is less than the specified 
        % tolerance.
        if abs(nCost_current - nCost_previous) < nEps
            break;
        end
        
        % Update cost. 
        nCost_previous = nCost_current;
    end

    % Return the optimized skew matrix.
    M_opt = M;

end % optimSkew


function dJ_dM = computePartialDerivative(Xdot, X, M)
%COMPUTEPARTIALDERIVATIVE Compute the partial derivative of a cost function 
% with respect to a matrix M using finite differences.
%
%   dJ_dM = computePartialDerivative(Xdot, X, M)
%
%   INPUT ===========================================================
%
%   Xdot (numeric matrix)
%   The rate of change of state variable X. 
%   Example: [1 0.5; 0.5 1]
%
%   X (numeric matrix)
%   The state variable matrix.
%   Example: [2 1; 1 2]
%
%   M (numeric matrix)
%   The matrix with respect to which the partial derivative will be computed.
%   Example: [0.5 1; 1 0.5]
%
%   OUTPUT ==========================================================
%
%   dJ_dM (numeric matrix)
%   The partial derivative of the cost function J with respect to matrix M.
%
%   PROCESS =========================================================
%
%   The function computes the partial derivative using the central difference
%   method. It perturbs each element of the matrix M by a small value epsilon
%   and evaluates the change in the cost function, J.
%
%   AUTHOR ==========================================================
%
%   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
%
%   =================================================================

% Set the step size for finite differences
epsilon = 1e-6;

% Get the size of Xdot to determine the dimensions of dJ_dM
[~, m] = size(Xdot);

% Initialize the partial derivative matrix
dJ_dM = zeros(m);

% Compute the partial derivative for each element of M
for i = 1:m
    for j = 1:m
        % Compute M_plus_eps by adding epsilon to the (i, j) element of M
        M_plus_eps        = M;
        M_plus_eps(i, j)  = M_plus_eps(i, j) + epsilon;

        % Compute M_minus_eps by subtracting epsilon from the (i, j) element of M
        M_minus_eps       = M;
        M_minus_eps(i, j) = M_minus_eps(i, j) - epsilon;

        % Compute J_plus_eps and J_minus_eps by evaluating the cost function
        J_plus_eps        = norm(Xdot - X*M_plus_eps, 'fro');
        J_minus_eps       = norm(Xdot - X*M_minus_eps, 'fro');

        % Compute the (i, j) element of the partial derivative matrix
        dJ_dM(i, j)       = (J_plus_eps - J_minus_eps) / (2 * epsilon);
    end
end
end % computePartialDerivative

%% ========================================================================
%        _       _   _   _                    
%  _ __ | | ___ | |_| |_(_)_ __   __ _        
% | '_ \| |/ _ \| __| __| | '_ \ / _` |       
% | |_) | | (_) | |_| |_| | | | | (_| |       
% | .__/|_|\___/ \__|\__|_|_| |_|\__, |       
% |_|                            |___/        
%   __                  _   _                 
%  / _|_   _ _ __   ___| |_(_) ___  _ __  ___ 
% | |_| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
% |  _| |_| | | | | (__| |_| | (_) | | | \__ \
% |_|  \__,_|_| |_|\___|\__|_|\___/|_| |_|___/



function visualize3DProjections(projectionData)
%VISUALIZE3DPROJECTIONS Visualize the projections on the first 3 jVec in 3D 
% space.
%
%   visualize3DProjections(projectionData)
%
%   INPUT =================================================================
%
%   projectionData (numeric matrix)
%   Data matrix where each column represents a component (jVec) and each row
%       is a data point.
%   Example: Random matrix (1000x3)
%
%   OUTPUT =================================================================
%
%   A 3D scatter plot visualization.
%
%   AUTHOR =================================================================
%
%   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
%
%   =======================================================================

    figure;
    scatter3(projectionData(:,1), projectionData(:,2), projectionData(:,3), 'filled');
    xlabel('1st jVec'); ylabel('2nd jVec'); zlabel('3rd jVec');
    title('3D Projections on First 3 jVec');
    grid on;

end % visualize3DProjections





function createCometPlots(projectionData)
%CREATECOMETPLOTS Create comet plots for various vector combinations.
%
%   createCometPlots(projectionData)
%
%   INPUT =================================================================
%
%   projectionData (numeric matrix)
%       Data matrix where each column represents a component (jVec) and each row
%       is a data point.
%       Example: Random matrix (1000x3)
%
%   OUTPUT =================================================================
%
%   A series of 2D comet plots showing trajectories between various jVec 
%   combinations.
%
%   AUTHOR =================================================================
%
%   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
%
%   =======================================================================

sFont = 'Arial';
nVecs = size(projectionData, 2);  % Determine the number of vectors
maxVal = max(abs(projectionData(:)));  % Maximum absolute value for axis scaling

% Calculate the velocity of data points for all pairs of vectors
velocity = sqrt(sum(diff(projectionData).^2,2));
normalized_velocity = (velocity - min(velocity)) / (max(velocity) - min(velocity));

% Creating a smooth-transition colormap from deep green to deep red
cmap = [linspace(0,1,256)', linspace(0.5,0,256)', linspace(0,0,256)'];

% Determine grid size for subplots based on the number of vector pairs
nPairs = nVecs * (nVecs-1) / 2;
gridRows = floor(sqrt(nPairs));  
gridCols = ceil(nPairs/gridRows);

% Initialize figure
figure; hold on; 

% Iterate through data points
for k = 2:size(projectionData, 1)
    pairIdx = 0;
    
    % Loop over all unique pairs of vectors
    for i = 1:nVecs
        for j = i+1:nVecs
            pairIdx = pairIdx + 1;
            
            % Choose subplot
            subplot(gridRows, gridCols, pairIdx); hold on;

            % If first iteration, set labels and properties
            if k == 2
                ax = gca;
                ax.XAxisLocation = 'origin';
                ax.YAxisLocation = 'origin';
                ax.XLim = [-maxVal maxVal];
                ax.YLim = [-maxVal maxVal];
                xlabel(['j',num2str(i)], 'FontName', sFont, 'FontSize', 11);
                ylabel(['j',num2str(j)], 'FontName', sFont, 'FontSize', 11);
                title(['j',num2str(i),' vs j',num2str(j)], 'FontName', sFont, 'FontSize', 11);
                ax.LineWidth = 1;
                ax.Color = 'none';
                ax.FontName = sFont;
                ax.FontSize = 11;
                grid off;
            end
            
            % Determine color and thickness for the current segment
            colorIdx = ceil(normalized_velocity(k) * 255) + 1;
            color = cmap(colorIdx, :);
            thickness = 1 + 3 * normalized_velocity(k);

            % Plot the current segment in the current subplot
            plot(projectionData(k-1:k,i), projectionData(k-1:k,j), ...
                'LineWidth', thickness, 'Color', color);
        end
    end
    
    % Pause to create an animation effect
    pause(0.01);
end

hold off;
end




function visualizePlanesProjections(projectionData)
%VISUALIZEPLANESPROJECTIONS Visualize all projections on all combinations 
%   of two jVec planes.
%
%   visualizePlanesProjections(projectionData)
%
%   INPUT =================================================================
%
%   projectionData (numeric matrix)
%   Data matrix where each column represents a component (jVec) and each row
%       is a data point.
%   Example: Random matrix (1000x3)
%
%   OUTPUT =================================================================
%
%   A series of 2D plots showing projections onto planes formed by various 
%       jVec combinations.
%
%   AUTHOR =================================================================
%
%   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
%
%   =======================================================================

    [~, numComponents] = size(projectionData);
    
    for i = 1:numComponents-1
        for j = i+1:numComponents
            figure;
            plot(projectionData(:,i), projectionData(:,j), 'o');
            xlabel(['jVec' num2str(i)]); ylabel(['jVec' num2str(j)]);
            title(['Projection on Plane formed by jVec' num2str(i) ' and jVec' num2str(j)]);
            grid on;
        end
    end
end
