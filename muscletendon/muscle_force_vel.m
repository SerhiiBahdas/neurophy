function nF = muscle_force_vel(nV, varargin)
%MUSCLE_FORCE_VEL Dimensionless muscle force-velocity relationship describes
% the dependence of the maximum isometric force output on the speed of muscle 
% contraction. 
%
%   nF = muscle_force_vel(nV)
%   nF = muscle_force_vel(nV, 'bPlot', 1)
%   nF = muscle_force_vel(nV, 'bPlot', 1, 'sMethod', 'Yakovenko')
%
%   INPUT =================================================================
%
%   nV (numeric array)
%   Velocity of the muscle contraction, [n.u.]. 
%   Example: -[0:0.05:1]; 
%
%   [OPTIONAL INPUT]
%
%   bPlot (boolean)
%   Visualize force-length relationship. 
%   Example: 1
%
%   sMethod (string)
%   The original Hill's approximation of the force-length relationship pro-
%   duces negative force at high contraction velocities, which is not physi-
%   ologically plausible [1]. We added an option to use Yakovenko's modifi-
%   cation of the force velocity relationship [2], suitable for modeling 
%   dynamics of the fast contraciton. 
%   Example: "Hill", "Yakovenko"
%   
%   OUTPUT ================================================================
%
%   nF (numeric array)
%   Muscle force, [n.u.].  
%
%   EXAMPLE ===============================================================
%
%   figure; 
%
%   nV = -[0:0.05:1]; 
%   nF_Hill = muscle_force_vel(nV, 'sMethod', 'Hill');
%   nF_Yakovenko = muscle_force_vel(nV, 'sMethod', 'Yakovenko');
%   
%   plot(nV, nF_Hill, nV, nF_Yakovenko); 
%   legend('Hill','Yakovenko')
%
%   title ('dimensionless muscle force-velocity relationship');
%   xlabel('contraction velocity, n.u.');
%   ylabel('force, n.u.')
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

%% SET DEFAULT VALUES. The calssical muscle force-velocity relationship de-
%  scribed in [1] is set as a default. 

% Default muscle force-velocity relationship approximation method. 
sMethod_default = 'Hill'; 

% By default, do not visualize force-length relationship. 
bPlot_default = 0; 

%% FETCH INPUTS. Fetch required contraction velocity and optional input. 

% Check if the input signal is numeric. 
if ~isnumeric(nV)
    error('Input vector must be numeric.')
end 

% Create an input parser object with default property values.
p = inputParser;

% Fetch the name of the method used to approximate force-length (optional).
addOptional(p,'sMethod',sMethod_default);

% Fetch the name of the method used to approximate force-length (optional).
addOptional(p,'bPlot',bPlot_default);

% Parse parameters. Assign them to a structure. 
parse(p,varargin{:}); p = p.Results; 

% Check if the directive to visualize data is boolean. 
if ~checkBoolean(p.bPlot)
    error('Directive to visualize plot must be boolean.');
end % checkBoolean

% For Hill's approximation. 
if p.sMethod == "Hill"

    % Maximal (normalized) isometric force.
    Fo = 1;

    % Additional coefficients. 
    a = 0.399; b = 0.331; 

    % Force-length relationship.
    fv = @(v) b*(Fo+a)./(-v+b)-0.4;

% For Yakovenko's approximation.
elseif p.sMethod == "Yakovenko"

    % Additional parameters. 
    a = 0.2;  b = -4.25; 

    % Force-length relationship.
    fv = @(v) (1-a*(v>0)).*(1-exp(b*v))./(1+exp(b*v))+1;

% If the method name was specified incorrectly display an error. 
else
    error('Check the spelling of the approximation method.')

end % p.sMethod

% Compute force. 
nF = fv(nV); 

% Plot tendon force-velocity relationship. 
if p.bPlot == true

    figure; 
    plot(nV, nF);
 
    title('dimensionless tendon force-velocity relationship');
    xlabel('contraction velocity, n.u.');
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