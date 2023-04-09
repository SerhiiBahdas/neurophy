function nSlopeList = slopearray(nSig, nRate, tWindow)
%SLOPEARRAY Computes linear trends in the cumsum of the input signal using
%sliding time window. 
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
%   tWindow (numeric)
%   Duration of the time window. 
%   Example: 10e-3
%
%   OUTPUT ============================================================
%   
%   nSlopeList (numeric array)
%   List of slopes computed in time windows.
%
%   AUTHOR ============================================================
%
%   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
%
%   ===================================================================


% Create time vector. 
tTime = signaltime(nSig, nRate); 

% Compute number of time windows signal can fit. 
numWindow = floor(tTime(end)/tWindow); 

% Preallocate memory. 
nSlopeList = zeros(1,numWindow); 

% Start sample. 
tWindow_start = 0; 

% End sample. 
tWindow_end = tWindow; 

    % Loop through windows. 
    for iWindow = 1:numWindow
    
        % Extract signal within the window. 
        nSig_window = extractSigPart(nSig, nRate, tWindow_start, tWindow_end); 
    
        % Compute slope of the linear trend in the signal. 
        nSlopeList(iWindow) = sigslope(nSig_window, nRate); 
    
        % Update start and end time of the window. 
        tWindow_start = tWindow_start + tWindow; 
        tWindow_end = tWindow_end + tWindow; 
    
    end % for

end % function