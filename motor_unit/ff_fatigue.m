function nForce = ff_fatigue(tTime, varargin)
%FF_FATIGUE Computes a force response of the fast fatiguable motor unit to
% repeated stimulation at a level that initially evokes maximum tension. 
%
%   nForce = ff_fatigue(tTime)
%   nForce = ff_fatigue(tTime, 'bPlot', 1)
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
%   Visualize computed force response of the fast fatiguable motor unit to
%   repeated stimulation at a level that initially evokes maximum tension. 
%   Example: 1
%
%   OUTPUT ================================================================
%
%   nForce (numeric array)
%   Maximal force output of the fatiguing motor unit, N.  
%
%   EXAMPLE ===============================================================
%
%   Simple:
%   tTime = [20:300]; 
%   nForce = ff_fatigue(tTime)
%
%   Advanced:
%   tTime = [10:500]; 
%   nForce = ff_fatigue(tTime, 'bPlot', 1)
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

% Compute first part of the muscle unit force (g) vs time (min) curve. 
% R-square = 0.9999

a1 = 176;
b1 = 0.3196;
c1 = 2.647;
a2 = 14.56;
b2 = 4.866;
c2 = -0.419;
a3 = 3.928;
b3 = 9.674;
c3 = 2.71;
a4 = 0.8212;
b4 = 14.74;
c4 = -0.8784;

nForce1 = a1*sin(b1.*tTime_min+c1) + a2*sin(b2.*tTime_min+c2) + ...
    a3*sin(b3.*tTime_min+c3) + a4*sin(b4.*tTime_min+c4); 

% Change unit of force to Newtons. 
nForce1 = nForce1/1000/9.8; 

% Compute second part of the muscle unit force (g) vs time (min) curve. 
% R-square = 0.9834.

a = 16.46;
b = -0.8629;
c = 2.164;
d = -0.04586;

nForce2 = a*exp(b*tTime_min) + c*exp(d*tTime_min); 

% Change unit of force to Newtons. 
nForce2 = nForce2/1000/9.8; 

% Find index of the sample where the first curve ends and the second one
% starts. 
[~,idx1] = min(abs(tTime_min-1.722095672)); 

% Cobine two curves into one. 
nForce             = zeros(1, length(nForce2)); 
nForce(1:idx1)     = nForce1(1:idx1); 
nForce(idx1+1:end) = nForce2(idx1+1:end); 

% Visualize data. 
if p.bPlot == 1

    figure; 

    % Connect computer points with a line. Add asterisks at datapoints.  
    plot(tTime, nForce, '-'); 

    % Add labels.
    xlabel('time, s'); 
    ylabel('force, N'); 
    title('Response of the FF MU to repeated stimulation at a level of that initially evokes maximum tension.')

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