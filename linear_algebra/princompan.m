function [nProj, nPC, nEVal, nVarExp] = princompan(nData, bPlot)
%PRINCOMPAN Principal component analysis on signal data.
%
%   [nProj, nPC, nEVal, nVarExp] = princompan(nData, bPlot)
%
%   Performs principal component analysis (PCA) on signal data for
%   dimensionality reduction and feature extraction. This function computes
%   the covariance matrix, performs eigenvalue decomposition, and calculates
%   the variance explained by each principal component.
%
%   INPUT =================================================================
%
%   nData (nSignals x nSamples numeric array):
%   Matrix of signal data where each row represents a signal and each
%   column represents a sample.
%
%   bPlot (boolean):
%   Enable plotting of the original signals, their principal components, and
%   the reconstruction of the original signals from the principal components.
%
%   OUTPUT ================================================================
%
%   nProj (numeric array):
%   Projections of the original data onto the principal components.
%
%   nPC (numeric array):
%   Principal components matrix where each column is a principal component.
%
%   nEVal (numeric array):
%   Sorted eigenvalues corresponding to the principal components.
%
%   nVarExp (numeric array):
%   Variance explained by each principal component.
%
%   EXAMPLE ===============================================================
%
%   % Generating Lorenz attractor data
%   sigma = 10; rho = 28; beta = 8/3;
%   lorenz = @(t, y)[sigma*(y(2)-y(1)); y(1)*(rho-y(3))-y(2); y(1)*y(2)-beta*y(3)];
%   [t, sol] = ode45(lorenz, [0 0.01:0.01:1], [1 1 1]);
%   [nProj, nPC, nEVal, nVarExp] = princompan(sol', true);
%
%   AUTHOR ================================================================
%
%   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
%
%   =======================================================================


% PCA ---------------------------------------------------------------------

% Computing covariance:
numSig           = size(nData,1);
nData_centered   = nData - mean(nData, 2);
nCov             = (nData_centered*nData_centered') / (numSig - 1);

% Eigenvalue decomposition:
[nEVec, nEVal]   = eig(nCov); 
[nEVal, idx]     = sort(diag(nEVal), 'descend');

% Principal components & projections: 
nPC              = nEVec(:, idx);
nProj            = nPC' * nData_centered;

% Calculate variance explained:
nTotalVar        = sum(nEVal);
nVarExp          = nEVal / nTotalVar;

% PLOTTING ----------------------------------------------------------------

% Plot if requested
if bPlot

    % Original Signals
    figure;
    subplot(3,1,1);
    plot((nData_centered + mean(nData, 2))');
    title('Original Signals');
    xlabel('Sample');
    ylabel('Amplitude');
    set(gca, 'LineWidth', 1.5, 'FontSize', 12, 'Box', 'off');
    
    % Projections
    subplot(3,1,2);
    plot(nProj');
    title('Projections');
    xlabel('Sample');
    ylabel('Amplitude');
    set(gca, 'LineWidth', 1.5, 'FontSize', 12, 'Box', 'off');
    
    % Reconstruction from PCs
    subplot(3,1,3);
    plot((nPC * nProj + mean(nData, 2))');
    title('Reconstruction from Principal Components');
    xlabel('Sample');
    ylabel('Amplitude');
    set(gca, 'LineWidth', 1.5, 'FontSize', 12, 'Box', 'off');

    % Enhance overall plot aesthetics
    set(gcf, 'Color', 'w');
    sgtitle('PCA Analysis and Reconstruction', 'FontSize', 14, ...
        'FontWeight', 'bold');
end

