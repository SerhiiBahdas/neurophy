function nSig_part = extractSigPart(nSig, nRate, tLimit_min, tLimit_max)
%EXTRACTSIGPART Extracts part of the signal from tLimit_min to tLimit_max. 
%
%   nSig_part = extractSigPart(nSig, nRate, tLimit_min, tLimit_max)
%
%   INPUT =============================================================
%   
%   nSig (numeric array)
%   Signal. 
%   Example: [1,2,3,4,5,6]
%
%   nRate (numeric)
%   Sampling rate, Hz. 
%   Example: [0, 1, 2]
%
%   tLimit_min (numeric)
%   Lower limit, s. 
%   Example: 1
%
%   tLimit_max (numeric)
%   Upper limit, s. 
%   Example: 2
%
%   OUTPUT ============================================================
%   
%   nSig_part (numeric)
%   Part of the signal.
%
%   AUTHOR ============================================================
%
%   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
%
%   ===================================================================

% Check if the lower limit is zero
if round(tLimit_min*nRate) == 0 
    nStart = 1; 
else
    nStart = round(tLimit_min*nRate); 

end

nSig_part = nSig(nStart:round(tLimit_max*nRate)); 

end