function [separatedData, jPCVecs] = joPCA(data, varargin)
%JOPCA Perform jPCA analysis. 
%
%   [separatedData, jPCVecs] = joPCA(data, varargin)
%
%   INPUT ===========================================================
%
%   data (structure)
%   A structure containing numbered fields (corresponding to the trial 
%   types) with data matrices. The matrices' size is nxm, where n is a
%   number of samples and m is the number of signals. 
%   Example: data(1).nX = rand(100,3); data(2).nX = rand(100,3);
%
%   numPC (numeric)
%   An even number of components to use in PCA and jPCA analyses.    
%   Example: 4
%
%   bPlanes (boolean)
%   Visualize all jPCA planes. 
%   Example: 0
%
%   b3D (boolean)
%   Project data onto three jPCA vectors. 
%   Example: 1
%
%   bAnimate (boolean)
%   Animate evolution of signals' projection onto the first two components.
%   Example: 1
%
%   bAttract (boolean)
%   If the point is outside the unit circle, attract it to the unit circle.
%   If it is inside the unit circle, repel it towards the radius of the 
%   unit circle. 
%   Example: 0
%
%   OUTPUT ==========================================================
%
%   separatedData (structure)
%   Input data projected onto the jPC space. 
%   
%   jPCVecs (numeric array)
%   Extracted jPC vectors. 
%
%   EXAMPLE =========================================================
%
%   % DATA MATRIX PARAMETERS.
%   numDim     = 10;    % number of dimensions
%   numSamples = 10000; % number of samples
% 
%   % SIGNAL PARAMETERS.
%   nA     = 1;     % semi-major axis
%   nB     = 0.5;   % semi-minor axis
%   nFreq  = 10;    % frequency
%   nPhase = pi/4;  % phase
% 
%   % GENERATE 2-D TRAJECTORY.
%   tTime = linspace(0, 5*pi, numSamples);
%   nTraj  = [nA * cos(nFreq * tTime); nB * sin(nFreq * tTime + nPhase)];
% 
%   % PROJECT SIGNALS TO 10-D.
%   nB = orth(randn(numDim, 2));
%   nX = (nB * nTraj)';
% 
%   % GENERATE DUMMY DATA FOR DIFFERENT 'TRIAL TYPES' OR 'CONDITIONS'.
%   dummydata(1).nX  = abs(1*nX);   dummydata(2).nX  = abs(2*nX); 
%   dummydata(3).nX  = abs(3*nX);   dummydata(4).nX  = abs(4*nX); 
%   dummydata(5).nX  = abs(5*nX);   dummydata(6).nX  = abs(6*nX); 
%   dummydata(7).nX  = abs(7*nX);   dummydata(8).nX  = abs(8*nX); 
%   dummydata(9).nX  = abs(9*nX);   dummydata(10).nX = abs(10*nX); 
%   dummydata(11).nX = abs(10*nX);  dummydata(12).nX = abs(11*nX); 
% 
%   % PARAMETERS FOR JPCA. 
%   numPC    = 2;  % number of components
%   bPlanes  = 1;  % visualize all jPC planes
%   b3D      = 1;  % visualize projections on 3 jPC vectors
%   bAnimate = 1;  % animate projections
%   bAttract = 0;  % attract point in jPCA space to unit circle 
% 
%   % COMPUTE & VISUALIZE PROJECTIONS. 
%   [separatedData, jPCVecs] = joPCA(dummydata, 'numPC', numPC,...
%       'bPlanes', bPlanes, 'bAnimate', bAnimate, 'b3D', b3D, ...
%       'bAttract', bAttract); 
%
%   AUTHOR ==========================================================
%
%   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
%
%   =================================================================

%% DEFAULT VALUES. Assign default values here. 

% Number of principal components and jPCs. 
numPC_default = 4; 

% Visualize planes formed by the jPCs. 
bPlanes_default = 0; 

% Animate evolution of signals' projection onto the first two components.
bAnimate_default = 0; 

% 3D plot option. 
b3D_default = 1;

% Attract point to a unit circle. 
bAttract_default = 0; 

%% FETCH ARGUMENTS. Fetch optional arguments and check their type. If argu-
% ment is not specified, assign a default value.  

% Create a structure to store inputs. 
p = inputParser; 

% Function that checks if the input is boolean. 
isboolean = @(x) islogical(x) || x == 1 || x == 0; 

% Add optional argument #1 to the structure.
addParameter(p,'numPC', numPC_default, @isnumeric);

% Add optional argument #2 to the structure.
addParameter(p, 'bPlanes', bPlanes_default, isboolean);

% Add optional argument #3 to the structure.
addParameter(p, 'bAnimate', bAnimate_default, isboolean);

% Add optional argument #4 to the structure.
addParameter(p, 'b3D', b3D_default, isboolean);

% Add optional argument #5 to the structure.
addParameter(p, 'bAttract', bAttract_default, isboolean);

% Assign the parameters to the structure.
parse(p,varargin{:}); p = p.Results; 

% Assign input arguments to variables to improve readability. 
numPC    = p.numPC;       bPlanes = p.bPlanes;      bAttract = p.bAttract; 
bAnimate = p.bAnimate;    b3D     = p.b3D;       

%% ORGANIZE DATA. Trials are concatenated together to capture their shared 
%  dynamics during the PCA and jPCA.

% Request even number of principal components. 
if mod(numPC,2)
    error('Specify even number of principal components (numPC).')
end % mod

% Get number of trial types. 
numTrialTypes = size(data,2); 

% Get data matrix name. 
sMatrixName = string(fieldnames(data(1))); 

% Initialize a variable to store the number of samples.
numSamples = 0;

% Iterate through the trial types in the structure.
for iTrialType = 1:numTrialTypes

    % Add the number of samples (rows) in the current data matrix to 
    % the variable. 
    numSamples = numSamples + size(data(iTrialType).(sMatrixName), 1);
    
end % iTrialType

% Get number of signals. 
numSigs = size(data(1).(sMatrixName),2); 

% Initialize an empty matrix to store the concatenated data.
X = zeros(numSamples, numSigs);

% Initialize an empty matrix to store start and end samples of each 
% trial type.
nStartEndSamples = zeros(numTrialTypes, 2);

% Create a counter.
iStartRow = 1;

% Loop through trial types.
for iTrialType = 1:numTrialTypes

    % Fetch data.
    nTmp = data(iTrialType).(sMatrixName);

    % Center data.
    nTmp = nTmp - mean(nTmp);

    % Compute number of rows.
    numRows = size(nTmp, 1);
   
    % Store the start and end samples of the current trial type.
    nStartEndSamples(iTrialType, 1) = iStartRow;
    nStartEndSamples(iTrialType, 2) = iStartRow + numRows - 1;

    % Assign data to the matrix.
    X(iStartRow:iStartRow + numRows - 1, :) = nan2zero(nTmp);

    % Increment counter.
    iStartRow = iStartRow + numRows;

end % iTrialType

% Assign nan-s to zero.
X(isnan(X)) = 0; 

%% JPCA. Find scores, skew-symmetric matrix, eigenvectors, and jPCs.

% Perform PCA on the input data
[~, PCs, ~] = pca(X, 'NumComponents', numPC, 'Centered', false, 'Algorithm', 'eig');

% Calculate derivative of the score matrix.
dPCs = diff(PCs, 1);

% Remove the last row from score matrix to match dimensions.
PCs = PCs(1:end-1, :); 

% Reflect trial type subtraction in the matrix storing sample #.
nStartEndSamples(end, 2) = nStartEndSamples(end, 2) - 1; 

% Calculate the least squares solution M.
M = pinv(PCs' * PCs) * PCs' * dPCs;

% Calculate skew-symmetric part of M.
Mskew = (M - M') / 2;

% Display skew-symmetric part of M.
disp('Skew-symmetric part of M:'); disp(Mskew);

% Optimize skew-symmetric transformation using gradient descent method. 
Mskew_opt = optimSkew(Mskew, PCs, dPCs);

% Display optimized skew-symmetric matrix.
disp('Optimized skew-symmetric matrix:'); disp(Mskew_opt);

% Perform eigen-decomposition of Mskew.
[V, D, ~] = eig(Mskew_opt);

% Get complex conjugate pairs of eigenvectors and eigenvalues.
[~, I] = sort(imag(diag(D)), 'descend'); V = V(:, I);

% Initialize jPC vectors.
jPCVecs = zeros(numPC);

% Calculate real-valued jPC vectors from complex eigenvectors.
for iComponent = 1:2:numPC
    jPCVecs(:, iComponent)   = real(V(:, iComponent)) - imag(V(:, iComponent+1));
    jPCVecs(:, iComponent+1) = real(V(:, iComponent)) + imag(V(:, iComponent+1));
end


%% PROJECTING. Project data from each trial type to jPC space. 

% Project the score matrix onto the jPCA plane.
nProj = PCs * jPCVecs;

% Initialize a structure to store data for each trial type data projection 
% separately.
separatedData = struct();

% Loop through trial types.
for iTrialType = 1:numTrialTypes
    
    % Get the start and end samples of the current trial type.
    iStartSample = nStartEndSamples(iTrialType, 1);
    iEndSample   = nStartEndSamples(iTrialType, 2);

    % Extract the segment corresponding to the current trial type.
    nSegment = nProj(iStartSample:iEndSample, :);
    
    % Assign the segment to the structure.
    separatedData(iTrialType).(sMatrixName) = nSegment;

end % iTrialType

%% VISUALIZE. Plot data projected on jPCA. 

% Minor bookkeeping for plotting projections onto jPC planes --------------

% Find all possible combinations between the two arrays.
[nC1, nC2] = meshgrid(1:numPC, 1:numPC); nCombs = [nC1(:) nC2(:)];

% Remove combinations with similar numbers.
nCombs = nCombs(nCombs(:,1) ~= nCombs(:,2), :);

% Remove duplicate combinations ([n,m] and [m,n] considered the same).
nCombs = unique(sort(nCombs, 2), 'rows');

% Number of remaining combinations. 
numCombs = length(nCombs); 

% -------------------------------------------------------------------------

% If the number of desired components is greater than two.  
if bPlanes && length(nCombs) > 2

    figure;

    % Visualize projections to all planes. 
    for iComp = 1:length(nCombs)

        % Loop through trial types. 
        for iTrialType = 1:numTrialTypes

            % Fetch the data for one trial type. 
            nData_onetrialtype = separatedData(iTrialType).(sMatrixName);

            % Subplot. 
            subplot(2, round(numCombs/2), iComp); hold on; 

            % Plot projections. 
            plot(nData_onetrialtype(:, nCombs(iComp,1)),...
                 nData_onetrialtype(:, nCombs(iComp,2)),...
                 'LineWidth', 0.2);

            % Add info.
            xlabel("jPC#" + string(nCombs(iComp,1)));
            ylabel("jPC#" + string(nCombs(iComp,2)));

        end % iTrialType
    end % iComb

    % Add a legend. 
    legend(); 

end % if

% If planes do not need to be visualized & user wants a static plot. 
if bPlanes==0 && bAnimate==0

    figure;

    % Loop through trial types. 
    for iTrialType = 1:numTrialTypes

        % Fetch the data for one trial type. 
        nData_onetrialtype = separatedData(iTrialType).(sMatrixName);

        % Plot projections. 
        plot(nData_onetrialtype(:, 1), nData_onetrialtype(:, 2),...
            'LineWidth', 0.2); hold on; 

        % Add info.
        title('Signals` projection onto the first two jPC vectors');
        xlabel("jPC#" + string(1));
        ylabel("jPC#" + string(2));

    end % iTrialType

    % Add a legend. 
    legend(); 

end %if

% If the number of desired components is greater than two and 3D option is selected.
if numPC > 2 && b3D
    figure;

    % Loop through trial types for PCA.
    for iTrialType = 1:numTrialTypes

        % Fetch the data for one trial type.
        nData_onetrialtype = PCs(nStartEndSamples(iTrialType, 1):nStartEndSamples(iTrialType, 2), :);

        % Subplot for PCA.
        subplot(1, 2, 1);
        plot3(nData_onetrialtype(:, 1), nData_onetrialtype(:, 2), nData_onetrialtype(:, 3),...
            'LineWidth', 0.2); hold on;

        % Add info.
        title('Signals` projection onto the first three PCA vectors');
        xlabel("PC#" + string(1));
        ylabel("PC#" + string(2));
        zlabel("PC#" + string(3));

    end % iTrialType

    % Loop through trial types for jPCA.
    for iTrialType = 1:numTrialTypes

        % Fetch the data for one trial type.
        nData_onetrialtype = separatedData(iTrialType).(sMatrixName);

        % Subplot for jPCA.
        subplot(1, 2, 2);
        plot3(nData_onetrialtype(:, 1), nData_onetrialtype(:, 2), nData_onetrialtype(:, 3),...
            'LineWidth', 0.2); hold on;

        % Add info.
        title('Signals` projection onto the first three jPC vectors');
        xlabel("jPC#" + string(1));
        ylabel("jPC#" + string(2));
        zlabel("jPC#" + string(3));

    end % iTrialType

    % Add a legend.
    legend();
end % if

% If user wants animated plot. 
if bAnimate==1

    % Visualize time evolution of each trial type's data projection onto
    % two jPC vectors. 
    cometPlot(separatedData, sMatrixName, bAttract); 
    
end % bAnimate

end % joPCA

%%                 _   _                      
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

            % Display every 10 iteration. 
            disp("iteration..." + string(iItter))

            % Reduce learning rate. 
            nLearningRate = nLearningRate/10; 

        end

        % If the
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
        M_plus_eps = M;
        M_plus_eps(i, j) = M_plus_eps(i, j) + epsilon;

        % Compute M_minus_eps by subtracting epsilon from the (i, j) element of M
        M_minus_eps = M;
        M_minus_eps(i, j) = M_minus_eps(i, j) - epsilon;

        % Compute J_plus_eps and J_minus_eps by evaluating the cost function
        J_plus_eps = norm(Xdot - X*M_plus_eps, 'fro');
        J_minus_eps = norm(Xdot - X*M_minus_eps, 'fro');

        % Compute the (i, j) element of the partial derivative matrix
        dJ_dM(i, j) = (J_plus_eps - J_minus_eps) / (2 * epsilon);
    end
end
end % computePartialDerivative

function A = nan2zero(A)
% NAN2ZERO Replace NaN values with zeros in an array
%
% Syntax: A = nan2zero(A)
%
% Inputs:
%   A - Input array
%
% Outputs:
%   A - Output array with NaN values replaced by zeros

A(isnan(A)) = 0;

end


%%       _       _   _   _                    
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

                                            
function cometPlot(separatedData, sMatrixName, bAttract)
% COMETPLOT Plot an evolution of the projections of data from each 
% trial type onto jPC vectors.
%
%   cometplot(sSeparatedData, sMatrixName)
%
%   INPUT ===========================================================   
%
%   sSeparatedData (structure)
%   A structure array containing data matrices. The structure must 
%       have at least 1 field (trial type).
%   Example: separatedData(1).nProximal = rand(100,5)
%
%   sMatrixName (string) 
%   Name of the matrix to plot.
%   Example: "nProximal"
%
%   AUTHOR ==========================================================
%
%   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
%
%   =================================================================

% Number of trial types (fields) in the structure.
numTrialTypes = numel(separatedData);

% Error is raised if the separatedData structure is empty.
if numTrialTypes < 1
    error('The input structure must contain at least 1 field.');
end

% Initialize variable storing max. number of samples from data fetched 
% from all trial types. 
nMaxSmpl = 0;

% Create a cell array to store matrices. 
matrices = cell(1, numTrialTypes);

% Initialize a vector containing information about the distance of the
% first sample from each trial type's data matrix from point [0,0]. 
nInitialDistances = zeros(1, numTrialTypes);

% Loop through trial types. 
for iTrialType = 1:numTrialTypes

    % Fetch matrices. 
    matrices{iTrialType} = separatedData(iTrialType).(sMatrixName);

    % Update a variable storing the maximum number of samples. 
    nMaxSmpl = max(nMaxSmpl, size(matrices{iTrialType}, 1));

    % Find distance of the initial value of each data matrix from the
    % origin [0,0] in the coordinate system formed by jPC1 and jPC2. 
    nInitialDistances(iTrialType) = norm(matrices{iTrialType}(1, 1:2));
end

% Normalize initial distances to the range [0, 1].
nNormalizedDistances = (nInitialDistances - min(nInitialDistances)) /...
    (max(nInitialDistances) - min(nInitialDistances));

% Assign nan-s to 0. 
nNormalizedDistances(isnan(nNormalizedDistances)) = 0; 

% Create figure and set up axes.
figure; ax = axes;
% ax.NextPlot = 'replaceChildren';

% Create a custom colormap. 
colormap(ax, customColors(numTrialTypes));

% Intialize axes with colormap. 
colors = colormap(ax);

% Length of the segment to visualize each iteration.
nSegment = round(nMaxSmpl/1000 + 1);

% Loop through data points. 
for iSample = 1:nSegment:nMaxSmpl

    % Loop through trials. 
    for iTrialType = 1:numTrialTypes

        % Fetch data matrix. 
        nMatrix = matrices{iTrialType};

        % Attract or repel points to the unit circle radius. 
        if bAttract == 1

            % Calculate the distance from the origin.
            nDist_origin = sqrt(nMatrix(:, 1).^2 + nMatrix(:, 2).^2);

            nMatrix(:, 1) = nMatrix(:, 1)./nDist_origin; 
            nMatrix(:, 2) = nMatrix(:, 2)./nDist_origin;       

        end % bAttract

        % If the current sample number does not exceed the maximum 
        % sample number.
        if iSample <= size(nMatrix, 1)

            % Find the ID of the color to use for the projection
            % visualization. 
            idColor = round(nNormalizedDistances(iTrialType) *...
                (size(colors, 1) - 1)) + 1;

            % Plot projection in the jPC vector space. 
            plot(nMatrix(1:iSample, 1), nMatrix(1:iSample, 2), ...
                'LineWidth', 2, 'Color', colors(idColor, :));
            hold on;

            % Add a marker of the spearhead of the projection. 
            plot(nMatrix(iSample, 1), nMatrix(iSample, 2), 'o',...
                'MarkerSize', 10, 'MarkerFaceColor',...
                colors(idColor, :), 'Color', colors(idColor, :));
        end % if iSample
    end % iTrialType

    % Unit circle radius. 
    if bAttract == 1
        xlim([-1,1])
        ylim([-1,1])

    else

        % Find min and max amplitude of the jPC1 and jPC2 observed in 
        % all trial types to set axes' limits. 
        xlim([min(cellfun(@(x) min(x(:, 1)), matrices)),...
             max(cellfun(@(x) max(x(:, 1)), matrices))]);

        ylim([min(cellfun(@(x) min(x(:, 2)), matrices)),...
           max(cellfun(@(x) max(x(:, 2)), matrices))]);

    end % bAttract

    % Add info. 
    xlabel('jPC1'); ylabel('jPC2');
    title('Evolution of the projections on the jPC vectors');
    drawnow; hold off;
    
    % Adjust the pause duration for faster or slower animation.
    pause(0.05); 

end % iSample
end % cometPlot


function nCmap = customColors(numColors)
%CUSTOMCOLORS Truncate standart colormap called 'hot' to remove white/
%bright yellow colors. Output n colors from the truncated colorscheme,
%where n is a number of trial types given as input. 
%
%   nCmap = customColors(numColors)
%
%   INPUT ===========================================================
%
%   numColors (numeric)
%   The number of colors to output. 
%   Example: 5
%
%   OUTPUT ==========================================================
%
%   nCmap (numeric)
%   Colormap with the desired number of colors.
%
%   AUTHOR ==========================================================
%
%   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
%
%   =================================================================

% Load standart colormap called 'hot'.
nCmap_hot = hot;

% Determine the number of colors in the original colormap.
numColor_orig = size(nCmap_hot, 1);

% Specify the fraction of colors you want to remove from the white
% and bright yellow end.
nFraction = 0.3; 

% Calculate the desired length of the new colormap.
numColor_new = round(numColor_orig * (1 - nFraction));

% Truncate the original colormap.
nCmap_truncated = nCmap_hot(1:numColor_new, :);

% Create the final colormap with the desired number of colors.
nCmap = interp1(1:numColor_new, nCmap_truncated,...
    linspace(1, numColor_new, numColors));

end % customColors