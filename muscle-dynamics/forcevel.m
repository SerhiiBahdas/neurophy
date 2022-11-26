function nF = forcevel(nV, varargin)
%FORCEVEL Force-velocity relationship describes the dependence of the
% maximum isometric force output on the speed of muscle contraction. 
%
%   nF = forcevel(nV)
%   nF = forcevel(nV, 'sMethod', 'Yakovenko')
%
%   INPUT =================================================================
%
%   nV (numeric array)
%   Velocity of the muscle contraction, m. 
%   Example: -[0:0.05:1]; 
%
%   [OPTIONAL INPUT]
%
%   sMethod (string)
%   The original Hill's approximation of the force-length relationship pro-
%   duces negative force at high contraction velocities that are not physi-
%   ologically plausible [1]. We added an option to use Yakovenko's modifi-
%   cation of the force velocity relationship [2], suitable for modeling 
%   dynamics of the fast contraciton. 
%   Example: "Hill", "Yakovenko"
%   
%   OUTPUT ================================================================
%
%   nF (numeric array)
%   Muscle force corresponding to the input contraction velocity.  
%
%   EXAMPLE ===============================================================
%
%   figure; 
%
%   nV = -[0:0.05:1]; 
%   nF_Hill = forcevel(nV, 'sMethod', 'Hill');
%   nF_Yakovenko = forcevel(nV, 'sMethod', 'Yakovenko');
%   
%   plot(nV, nF_Hill, nV, nF_Yakovenko); 
%   legend('Hill','Yakovenko')
%
%   title ('normalized muscle force-length relationship');
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

%% FETCH INPUTS. Fetch required contraction velocity and optional input. 

% Check if the input signal is numeric. 
if ~isnumeric(nV)
    error('Input vector is not numeric.')
end 

% Create an input parser object with default property values.
p = inputParser;

% Fetch the name of the method used to approximate force-length (optional).
addOptional(p,'sMethod',sMethod_default);

% Parse parameters. Assign them to a structure. 
parse(p,varargin{:}); p = p.Results; 

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

end % function