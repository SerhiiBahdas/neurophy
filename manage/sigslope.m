function nSlope = sigslope(nNoise, nRate)
%SIGSLOPE Computes slope of the linear trend in cummulative sum of the 
% signal.  
%
%   nSlope = noiseslope(nNoise, nRate)
%
%   INPUT =============================================================
%   
%   nNoise (numeric array)
%   Noise. 
%   Example: [1,2,3,4,5,6]
%
%   nRate (numeric)
%   Sampling rate. 
%   Example: [0, 1, 2]
%
%   OUTPUT ============================================================
%   
%   nSlope (numeric)
%   Slope of the linear trend.
%
%   AUTHOR ============================================================
%
%   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
%
%   ===================================================================

% Compute time vector. 
tTime = signaltime(nNoise, nRate); 

% Find linear trend in the signal (y ~ mx).
nSlope = polyfit(tTime, cumsum(nNoise), 1); 

% Save the slope of the linear trend. 
nSlope = nSlope(1); 

end