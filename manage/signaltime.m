function tTime = signaltime(nSig, nRate)
%SIGNALTIME Computes time vector for the intput signal, given its sampling
%frequency. 
%
%   tTime = signaltime(nSig, nRate)
%
%   INPUT =============================================================
%   
%   nSig (numeric array)
%   Input signal. 
%   Example: [1,2,3,4,5,6]
%
%   nRate (numeric)
%   Sampling rate. 
%   Example: [0, 1, 2]
%
%   OUTPUT ============================================================
%   
%   tTime (numeric array)
%   Time vector.
%
%   AUTHOR ============================================================
%
%   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
%
%   ===================================================================

tTime = 0:1/nRate:(length(nSig)-1)/nRate; 

end