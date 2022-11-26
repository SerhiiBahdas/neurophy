function nSig = emgfilt(nSig, nRate, varargin)
%EMGFILT Filters, demeans, rectifies, normalizes EMG, and computes its 
%   linear envelope. 
%
%   EXAMPLE ===============================================================
%
%   [STANDART PROCESSING]:
%
%   nSig = emgfilt(nSig, nRate)
%
%   [CUSTOM PROCESSING]:
%
%   nSig = emgfilt(nSig, nRate, 'nFreq_BPC', [20, 500], 'nOrder', 4,...
%                  'nFreq_LPC', 10, 'nMVC', 0.1)
%
%   INPUT =================================================================
%
%   nSig (numeric array)
%   Raw EMG signal. 
%   Example: [1,1,1,1,1]
%
%   nRate (double)
%   Sampling rate used to record EMG, Hz. 
%   Example: 2e3
%
%   [OPTIONAL INPUTS]
%
%   nFreq_BPC (numeric array)
%   Cutoff frequencies for Butterword bandpass filter, Hz.
%   Example: [20, 500]
%
%   nMVC (double)
%   Voltage value (many use maximal voluntary contraction) used to 
%       normalize EMG, V.
%   Example: 0.1
%
%   nFreq_LPC (double)
%   Cutoff frequency for Butterword lowpass filter used to create a 
%       linear envelope from normalized EMG. 
%   Example: 10
%
%   nOrder (double)
%   Order of the bandpass and lowpass Butterword filters. 
%   Example: 4
%
%   OUTPUT ================================================================
%
%   nSig (numeric array)
%   Linear envelope (bandpass filtred, demeaned, rectified, normalized, 
%       and lowpass filtered EMG). 
%
%   AUTHOR ================================================================
%
%   Serhii Bahdasariants, WVU, NEL, https://github.com/SerhiiBahdas
%
%   REFERENCES ============================================================
%
%   1. Robertson, D. Gordon E., Graham E. Caldwell, Joseph Hamill, Gary 
%       Kamen, and Saunders N. Whittlesey. 2014. Research Methods in 
%       Biomechanics. Second edition. Champaign, Illinois: Human Kinetics.
%
%   =======================================================================

%% SET DEFAULT VALUES. The cutoff values used to low- and highpass filter 
%  the EMG are fetched from [1]. The default maximal voluntary contraction
%  value is calculated as the maximum amplitude of the input EMG signal. 

% Default cutoff frequencies for a bandpass filter, Hz. 
nFreq_BPC_default = [20, 500]; 

% Maximal voluntary activation (used to normalize EMG), V.
nMVC_default = max(nSig); 

% Default low-pass filter cutoff frequency (Hz) to create linear envelope
% from rectified EMG. 
nFreq_LPC_default = 10; 

% Default value of the Butterworth filter's order. 
nOrder_default = 4; 

%% FETCH INPUTS. Fetch required input EMG signal [1*n] and the sampling 
%  rate to record the EMG. Optional inputs include band- and lowpass
%  filters cutoff frequencies, order of the filters, and maximal voluntary 
%  activation value (used to normalize EMG). 

% Check if the input signal is numeric. 
if isnumeric(nSig)

    % Convert to double (function filtfilt used later in the program 
    % requires input signals to be double). 
    nSig = double(nSig); 

else 
    % Otherwise, raise error. 
    error('Input signal is not numeric.')
end 

% Check if the input sampling rate is numeric. 
if ~isnumeric(nRate)
    error('Input sampling rate is not numeric.')
end 
    
% Create an input parser object with default property values.
p = inputParser;

% Fetch the cutoff frequencies for bandpass filter (optional).
addOptional(p,'nFreq_BPC',nFreq_BPC_default,@isnumeric);

% Fetch the cutoff frequency for lowpass filter used to create linear
% envelope from rectified EMG (optional). 
addOptional(p,'nFreq_LPC',nFreq_LPC_default,@isnumeric);

% Fetch order of the butterworth filter (optional). 
addOptional(p,'nOrder',nOrder_default,@isnumeric);

% Fetch a value (many use maximal voluntary contraction) used to normalize 
% EMG, V (optional).
addOptional(p,'nMVC',nMVC_default,@isnumeric);

% Parse parameters. Assign them to a structure. 
parse(p,varargin{:}); p = p.Results; 

%% PROCESS EMG. Bandpass filter the signal to remove artifacts. Demean and 
%  normalize it. Lowpassfilter the resulting signal to create linear enve-
%  lope. 

% Return the transfer function coefficients of an nth-order bandpass 
% digital Butterworth filter with normalized cutoff frequency.
[b1,a1] = butter(p.nOrder, p.nFreq_BPC/(nRate*0.5), 'bandpass'); 

% Use the bandpass filter to remove high frequency noise and prevent 
% aliasing (an effect that can make different signals indistiguishable 
% from one another), and baseline drift (movement artifacts, moisture on 
% the skin, and DC offset).
nSig = filtfilt(b1,a1,nSig); 

% Demean bandpass-filtered EMG signal. 
nSig = nSig - mean(nSig); 

% Rectify the demeaned signal. 
nSig = abs(nSig); 

% Normalize EMG amplitude. 
nSig = nSig./p.nMVC; 

% Return the transfer function coefficients of an nth-order lowpass digital
% Butterworth filter with normalized cutoff frequency.
[b2,a2] = butter(p.nOrder, p.nFreq_LPC/(nRate*0.5), 'low'); 

% Lowpass filter the normalized signal to create linear envelope. 
nSig = filtfilt(b2,a2,nSig); 

end % emgfilt
