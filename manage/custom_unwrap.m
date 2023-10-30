function unwrapped = custom_unwrap(phase, tolerance)
%CUSTOM_UNWRAP Unwrap phase data considering a given tolerance around pi.
%
%   UNWRAPPED = CUSTOM_UNWRAP(PHASE, TOLERANCE) unwraps the phase data by adjusting
%   phase jumps that are close to pi, based on a specified tolerance.
%
%   INPUT =================================================================
%
%   PHASE (numeric array)
%       Input phase data. An array containing phase values.
%
%   TOLERANCE (numeric scalar)
%       Tolerance around pi for unwrapping. If the phase difference between 
%       consecutive elements is within [pi - TOLERANCE, pi + TOLERANCE], 
%       the phase will be unwrapped.
%
%   OUTPUT ================================================================
%
%   UNWRAPPED (numeric array)
%       Unwrapped phase data. An array of phase values that have been adjusted
%       to account for phase jumps.
%
%   FUNCTION OPERATIONS ===================================================
%
%   (1) Initializes the unwrapped phase with the first value of input phase.
%   (2) Iterates through the phase data, comparing consecutive elements.
%   (3) Adjusts the current phase value based on the difference from the previous
%       value and the specified tolerance.
%
%   AUTHOR ================================================================
%
%   ChatGPT, OpenAI, https://openai.com
%
%   =======================================================================

    if nargin < 2
        tolerance = 0.1; % Default tolerance around pi
    end

    % Ensure the phase data is a column vector
    phase = phase(:);

    % Initialize the unwrapped phase with the first value of input phase
    unwrapped = zeros(size(phase));
    unwrapped(1) = phase(1);

    for ii = 2:length(phase)
        delta = phase(ii) - phase(ii-1);
        
        while delta > pi + tolerance
            phase(ii) = phase(ii) - 2*pi;
            delta = phase(ii) - phase(ii-1);
        end

        while delta < -(pi + tolerance)
            phase(ii) = phase(ii) + 2*pi;
            delta = phase(ii) - phase(ii-1);
        end
        
        unwrapped(ii) = phase(ii);
    end
end
