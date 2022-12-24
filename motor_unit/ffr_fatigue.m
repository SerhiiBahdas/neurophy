function nForce = ffr_fatigue(tTime, varargin)
%FFR_FATIGUE Computes a force response of the fast fatigue-resistant motor 
% unit to repeated stimulation at a level that initially evokes maximum 
% tension. 
%
%   nForce = ffr_fatigue(tTime)
%   nForce = ffr_fatigue(tTime, 'bPlot', 1)
%
%   INPUT =================================================================
%
%   tTime (numeric)
%   Time elased after the onset of maximal contraction, s. 
%   Example: [20:300]; 
%
%   [OPTIONAL INPUT]
%
%   bPlot (boolean)
%   Visualize computed force response of the fast fatigue-resistant motor 
%   unit to repeated stimulation at a level that initially evokes maximum 
%   tension. 
%   Example: 1
%
%   OUTPUT ================================================================
%
%   nForce (numeric array)
%   Maximal force output of the fatiguing motot unit, N.  
%
%   EXAMPLE ===============================================================
%
%   Simple:
%   tTime = [20:300]; 
%   nForce = ffr_fatigue(tTime)
%
%   Advanced:
%   tTime = [10:500]; 
%   nForce = ffr_fatigue(tTime, 'bPlot', 1)
% 
%   REFERENCES ============================================================
%
%      1. Burke, R. E., D. N. Levine, P. Tsairis, and F. E. Zajac. 1973. 
%      "Physiological Types and Histochemical Profiles in Motor Units 
%      of the Cat Gastrocnemius." The Journal of Physiology 234 (3): 723
%      –48. https://doi.org/10.1113/jphysiol.1973.sp010369.
%
%   AUTHOR ================================================================
%   
%   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
%
%   =======================================================================

%% SET DEFAULT VALUES. By default, do not visualize the tension in reponse 
% to repretitive stimulation. 

% By default, do not visualize force-time relationship. 
bPlot_default = 0; 

%% FETCH INPUTS. Create input parser and fetch inputs. 

% Create an input parser object with default property values.
p = inputParser;

% Fetch the option to visualize plot. 
addOptional(p,'bPlot',bPlot_default);

% Parse parameters. Assign them to a structure. 
parse(p,varargin{:}); p = p.Results; 

% Check if the directive to visualize data is boolean. 
if ~checkBoolean(p.bPlot)
    error('Directive to visualize plot must be boolean');
end % checkBoolean

% Check if the input time vector is numeric. 
if ~isnumeric(tTime)
    error('Input timevector must be numeric.')
end 

% Measure time in minutes. 
tTime_min = tTime/60; 

% Compute muscle unit force (N) at a given momemnt of time (s). 
% R-square = 0.994.
nForce1 = -0.03243.*tTime_min.^5 +...
            0.4987.*tTime_min.^4 +...
           -2.94.*tTime_min.^3 +...
            7.817.*tTime_min.^2 +...
           -9.257.*tTime_min    +...
            34.39; 

% Change unit of force to Newtons. 
nForce1 = nForce1/1000/9.8; 

% Compute muscle unit force (g) at a given momemnt of time (m). 
% R-square = 0.994.
nForce2 = 49.23*exp(-0.2444.*tTime_min) + 13.79*exp(-0.02702.*tTime_min); 

% Change unit of force to Newtons. 
nForce2 = nForce2/1000/9.8; 

% Find index of the sample where the first curve ends and the second one
% starts. 
[~,idx1] = min(abs(tTime_min-5)); 

% Cobine two curves into one. 
nForce             = zeros(1, length(nForce2)); 
nForce(1:idx1)     = nForce1(1:idx1); 
nForce(idx1+1:end) = nForce2(idx1+1:end); 

% Visualize data. 
if p.bPlot == 1

    figure
    % Connect computer points with a line. 
    plot(tTime, nForce, '-'); 

    % Add labels.
    xlabel('time, s'); 
    ylabel('force, N'); 
    title('Response of the FFR MU to repeated stimulation at a level of that initially evokes maximum tension.')

end % p.bPlot

end % ff_fatigue

%% ADDITIONAL FUNCTIONS. Helper functions used in the script are listed below. 

function bBoolean = checkBoolean(bInput)
% CHECKBOOLEAN Checks if the input variable is boolean. 

    if (bInput == true) || (bInput == false)
        bBoolean = 1; 
    else
        bBoolean = 0;
    end % bInput

end % checkBoolean
