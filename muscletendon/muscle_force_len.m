function nF = muscle_force_len(nX, varargin)
%MUSCLE_FORCE_LEN Dimensionless muscle force-length relationship describes 
% dependence of the produced muscle force on its length. 
%
%   nF = muscle_force_len(nX)
%   nF = muscle_force_len(nX, 'bPlot', 1)
%   nF = muscle_force_len(nX, 'bPlot', 1, 'nOptLen', 0.5,...
%                         'sMovementType', 'dynamic')
%
%   INPUT =================================================================
%
%   nX (numeric array)
%   Muscle length, [n.u.]. 
%   Example: [0:0.05:1]; 
%
%   [OPTIONAL INPUT]
%
%   bPlot (boolean)
%   Visualize force-length relationship.
%   Example: 0
%
%   nOptLen (numeric)
%   Optimal muscle fiber length, [n.u.].
%   Example: 0.5
%
%   sMovementType(string)
%   Movement type. 
%   Study [1] showed that classical isometric length-tension curves of 
%   active wrist muscles [2] are not representative of continuous (dynamic)
%   movements. Thus, we added option to model force-length curve with no 
%   descending limb, as described in [3]. 
%   Example: 'static', 'dynamic'
%   
%   OUTPUT ================================================================
%
%   nF (numeric array)
%   Muscle force, [n.u.].  
%
%   EXAMPLE ===============================================================
%
%   nX = [0:0.05:1]; 
%   nF_static  = muscle_force_len(nX, 'bPLot', 1);
%   title ('static dimensionless muscle force-length relationship');    
%
%   nF_dynamic = muscle_force_len(nX, 'bPlot', 1, 'sMovementType', 'dynamic');
%   title ('dynamic dimensionless muscle force-length relationship');  
% 
%   REFERENCES ============================================================
%
%   1. Gillard, Deborah M, Sergiy Yakovenko, Tracy Cameron, and Arthur 
%      Prochazka. 2000. "Isometric Muscle Length–Tension Curves Do Not 
%      Predict Angle–Torque Curves of Human Wrist in Continuous Active 
%      Movements." Journal of Biomechanics 33 (11): 1341–48. 
%
%   2. Zajac, F. E. 1989. "Muscle and Tendon: Properties, Models, Scaling, 
%      and Application to Biomechanics and Motor Control." Critical Reviews
%      in Biomedical Engineering 17 (4): 359–411.
%
%   3. Yakovenko, S., V. Gritsenko, and A. Prochazka. 2004. "Contribution 
%      of Stretch Reflexes to Locomotor Control: A Modeling Study.” 
%      Biological Cybernetics 90 (2): 146–55. 
%
%   4. Thelen, Darryl G. 2003. "Adjustment of Muscle Mechanics Model 
%      Parameters to Simulate Dynamic Contractions in Older Adults.” 
%      Journal of Biomechanical Engineering 125 (1): 70–77.
%
%   5. Romero, F., and F. J. Alonso. 2016. “A Comparison among Different 
%      Hill-Type Contraction Dynamics Formulations for Muscle Force 
%      Estimation." Mechanical Sciences 7 (1): 19–29. 
%
%   AUTHOR ================================================================
%   
%   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
%
%   =======================================================================

%% SET DEFAULT VALUES. The calssical muscle force-length relationship with
%  the descending limb (valid for static postures) is set as default. Opt
%  the default optimal muscle length to be equal to 0.5 n.u. 

% Choose default method for approximating force-length realtionship. 
sMovementType_default = 'static'; 

% Optimal length of the muscle fiber, [n.u.]. 
nOptLen_default = 0.5; 

% By default, do not visualize force-length relationship. 
bPlot_default = 0; 

%% FETCH INPUTS. Fetch required contraction velocity and optional input. 

% Check if the input signal is numeric. 
if ~isnumeric(nX)
    error('Input vector of muscle lengths must be numeric.')
end 

% Create an input parser object with default property values.
p = inputParser;

% Fetch directive to visualize force-length relationship.
addOptional(p,'nOptLen',nOptLen_default);

% Fetch the type of the movement to know which approximation to use.
addOptional(p,'sMovementType',sMovementType_default);

% Fetch optimal length of the muscle fiber, [n.u.].
addOptional(p,'bPlot',bPlot_default);

% Parse parameters. Assign them to a structure. 
parse(p,varargin{:}); p = p.Results; 

% Check if the muscle fiber length is numerical value. 
if ~isnumeric(p.nOptLen)
    error('Input optimal muscle fiber length must be numeric.')
end 

% Check if the directive to visualize data is boolean. 
if ~checkBoolean(p.bPlot)
    error('Directive to visualize plot must be boolean.');
end % checkBoolean

%% COMPUTE FORCE. For a static condition use default function. For dynamic 
% (continuous) movements, use [3]. 

% Describe passive part of the contractile component muscle force-length
% relationship [4].
fl_passive = @(x)(exp(2*(x-p.nOptLen))-1)/(exp(1)-1).*(x>=p.nOptLen); 

% For static condition.
if p.sMovementType == "static"

    % Describe active part of the contractile component muscle force-length 
    % relationship [4,5]. 
    fl_active = @(x) exp((-(x - p.nOptLen).^2)/0.45); 

% For dynamic condition.
elseif p.sMovementType == "dynamic"

    % Describe active part of the contractile component muscle force-length 
    % relationship [1,3].
    fl_active = @(x) 2*x.*(x<p.nOptLen) + ((x-p.nOptLen)+1).*(x>=p.nOptLen);

% If the movement name was specified incorrectly display an error. 
else
    error('Check the spelling of the movement type; options: static, dynamic');

end % p.sMovementType

% Compute force proportional to muscle length. 
nF = fl_active(nX) + fl_passive(nX); 

%% VISUALIZE PLOT. Plot active, passive, and resulting force.

% Plot figure. 
if p.bPlot == true

    figure; 
    plot(nX, fl_active(nX), nX, fl_passive(nX), nX, nF); 
    
    legend('active component', 'passive component', 'resulting force'); 
    title ('dimensionless muscle force-length relationship');
    xlabel('muscle length, n.u.');
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