function [nMean, nSTD] = plotWstats(nSig, sColor)
%PLOTWSTATS Plot (1) signals from the input 2-D matrix, (2) average, (3)
%positive and negative standart deviations. The area between the standart
%deviations is shaded automatically. 
%
%   plotWstats(nSig, sColor)
%
%   INPUT =================================================================
%
%   nSig (numeric 2-D array)
%   Input matrix. 
%   Example: rand(100,3)
%
%   sColor (string)
%   Color. 
%   Example: 'k'
%
%   OUTPUT ================================================================
%   
%   nMean (numeric array)
%   Average of the input signals. 
%
%   nSTD (numeric array)
%   Standart deviation of the input signals. 
%
%   AUTHOR ================================================================
%   
%   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
%
%   =======================================================================

% Find the greater dimention. 
[row, col] = size(nSig); 

% Signals are in columns. Samples are in rows. 
if col > row

    % Transpose the matrix. 
    nSig = nSig'; 

end

% Compute average. 
nMean = mean(nSig, 2);

% Compute STD. 
nSTD = std(nSig, 0, 2); 

% Plot positive STD. 
nSTD_upc = nMean + nSTD; 
plot(nSTD_upc, 'LineWidth', 2, 'Color', sColor);

% Keep plotting.
hold on; 

% Plot negative STD. 
nSTD_dwc = nMean - nSTD; 
plot(nSTD_dwc, 'LineWidth', 2, 'Color', sColor);

% Shade the area between + and - STD. 
nXaxis = [1:col, flip(1:col)]; 
nYaxis = [nSTD_dwc', fliplr(nSTD_upc')]; 
fill(nXaxis', nYaxis', sColor);

% Plot average. 
plot(nMean, 'LineWidth', 3, 'Color', 'k'); 

% Plot signals. 
plot(nSig, 'LineWidth', 0.2, 'Color', sColor);

% Add transparency. 
alpha(.5); 

end % function