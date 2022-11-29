function nF = tendon_force_len(nX, varargin)
%TENDON_FORCE_LEN Dimensionless tendon force-length relationship describes 
% force output of a tendon as a function of its length based on its material 
% properties. 
%
%   nF = tendon_force_len(nX)
%   nF = tendon_force_len(nX, 'bPlot', 1)
%
%   INPUT =================================================================
%
%   nX (numeric array)
%   Tendon length, [n.u.]. 
%   Example: [0:.001:0.05]; 
%
%   [OPTIONAL INPUT]
%   
%   bPlot (boolean)
%   Visualize force-length relationship. 
%   Example: 1
%
%   OUTPUT ================================================================
%
%   nF (numeric array)
%   Tendon force, [n.u.].  
%
%   EXAMPLE ===============================================================
%
%   nX = [0:.001:0.05]; 
%   nF = tendon_force_len(nX, 'bPlot', 1);
% 
%   REFERENCES ============================================================
%
%   1. Zajac, F. E. 1989. "Muscle and Tendon: Properties, Models, Scaling, 
%      and Application to Biomechanics and Motor Control." Critical Reviews
%      in Biomedical Engineering 17 (4): 359–411.
%
%   2. Yakovenko, S., V. Gritsenko, and A. Prochazka. 2004. “Contribution 
%      of Stretch Reflexes to Locomotor Control: A Modeling Study.” 
%      Biological Cybernetics 90 (2): 146–55. 
%
%   AUTHOR ================================================================
%   
%   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
%
%   =======================================================================

%% SET DEFAULT VALUES. By default, do not visualize tendon force-length 
% relationship. 

% By default, do not visualize force-length relationship. 
bPlot_default = 0; 

%% FETCH INPUTS. Create input parser and fetch inputs. 

% Create an input parser object with default property values.
p = inputParser;

% Fetch the name of the method used to approximate force-length (optional).
addOptional(p,'bPlot',bPlot_default);

% Parse parameters. Assign them to a structure. 
parse(p,varargin{:}); p = p.Results; 

% Check if the directive to visualize data is boolean. 
if ~checkBoolean(p.bPlot)
    error('Directive to visualize plot must be boolean');
end % checkBoolean

% Check if the muscle fiber contraction velocity is numerical value. 
if ~isnumeric(nX)
    error('Input tendon length must be numeric.')
end 

%% COMPUTE FORCE. Compute tendon force with respect to its length. 

% According to [1,2], the tendon model can be described as follows.
fl_ten = @(x) 0.10377*(exp(91*x)-1).*(x<=0.01516) + ...
    (37.526*x-0.26029).*(x>0.01516);

% Compute output tendon force. 
nF = fl_ten(nX); 

% Plot tendon force-length relationship. 
if p.bPlot == true

    figure; 
    plot(nX, nF);
 
    title ('dimensionless tendon force-length relationship');
    xlabel('tendon length, n.u.');
    ylabel('force, n.u.')

end % p.bPlot 

end % function

%% ADDITIONAL FUNCTIONS. Add a helper function that checks if the variable 
%  type is boolean. This function is needed to catch errors in the main 
%  function. 

function bBoolean = checkBoolean(bInput)

    if (bInput == true) || (bInput == false)
        bBoolean = 1; 
    else
        bBoolean = 0;
    end % bInput

end % checkBoolean