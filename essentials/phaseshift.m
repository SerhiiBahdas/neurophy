function nSig = phaseshift(nSig, nTol)
%PHASESHIFT Unwraps the radian phase angles in a vector nSig. Whenever the
% jump between consecutive angles is greater than or equal to -π, π, -π/2, 
% or π/2 radians, PHASESHIFT shifts the angles by adding multiples of ±π or
% ±π/2.
%
%   nSig = phaseshift(nSig, nTol)
%
%   INPUT =================================================================
%
%   nSig (numeric array)
%   Input signal. 
%   Example: sin(1:10)
%
%   nTol (double)
%   Tollerance.
%   Example: 1e-2
%   
%   OUTPUT =========================================================
%
%   nSig (numeric array)
%   Processed signal.
%
%   AUTHOR =========================================================
%
%   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
%
%   ================================================================

% Number of samples
numSample = numel(nSig); 

% For all samples
for iSample = 2:numSample
    % Check if the offset is equal to π
    if approxequal(nSig(iSample) - nSig(iSample-1), pi, nTol)
        nSig(iSample:end) = nSig(iSample:end) - pi; 
    % Check if the offset is equal to -π
    elseif approxequal(nSig(iSample) - nSig(iSample-1), -pi, nTol)
        nSig(iSample:end) = nSig(iSample:end) + pi; 
    % Check if the offset is equal to π/2
    elseif approxequal(nSig(iSample) - nSig(iSample-1), pi/2, nTol)
        nSig(iSample:end) = nSig(iSample:end) - pi/2; 
    % Check if the offset is equal to -π/2
    elseif approxequal(nSig(iSample) - nSig(iSample-1), -pi/2, nTol)
        nSig(iSample:end) = nSig(iSample:end) + pi/2; 
    end % if

end % function