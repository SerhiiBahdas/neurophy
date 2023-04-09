function nSig_denoised = wavefilt(nSig, nRate, bPlot)
% WAVEFILT denoise EMG signal. 
%
%   EXAMPLE ===============================================================
%
%   nRate = 1000; 
%   nSig  = sin(0:1/nRate:10); 
%   nSig = nSig + randn(1,length(nSig))/20; 
%   bPlot = 1; 
%
%   nSig_denoised = wavefilt(nSig, nRate, bPlot); 
%
%   INPUT =================================================================
%
%   nSig (numeric array)
%   Raw EMG signal. 
%   Example: [1,1,1,1,1]
%
%   nRate (double)
%   Sampling rate, Hz.
%   Example: 2e3
%
%   bPlot (boolean)
%   Enable plotting. 
%   Example: 0
%
%   OUTPUT ================================================================
%
%   nSig_denoised (numeric array)
%   Denoised EMG. 
%
%   AUTHOR ================================================================
%
%   Serhii Bahdasariants, WVU, NEL, https://github.com/SerhiiBahdas
%
%   =======================================================================

% Minimal energy contribution of the component, %. 
nContrib = 0.05; 

% Name of the wavelet function. 
sWaveletName = 'sym8';

% Maximal overlap discrete wavelet packet transform. 
[wpt,~,~,~,nRelenergy] = modwpt(nSig, sWaveletName); 

% Levels, the relative energy of which is > nContrib. 
bLevelForReconstruction = nRelenergy > nContrib;

% Sum down the rows of the selected signals.
nSig_denoised = sum(wpt(bLevelForReconstruction,:),1);

% When plotting is enabled.
if bPlot == 1

    % Open figure. 
    figure;

    % Plot noisy and denoised signals. 
    plot(0:1/nRate:(length(nSig)-1)/nRate, nSig, 'b',...
         0:1/nRate:(length(nSig_denoised)-1)/nRate, nSig_denoised, 'r',...
         'LineWidth', 1); 

    % Add legend and labels. 
    legend(["original", "denoised"]); xlabel('time,s'); ylabel("EMG, V"); 

end % function












