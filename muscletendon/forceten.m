function nF = forceten(nX)
%FORCETEN describes force output of a tendon as a function of its length
%based on its material properties. 
%
%   nF = forceten(nX)
%
%   INPUT =================================================================
%
%   nX (numeric array)
%   Tendon length, [n.u.]. 
%   Example: [0:.001:0.05]; 
%
%   OUTPUT ================================================================
%
%   nF (numeric array)
%   Tendon force, [n.u.].  
%
%   EXAMPLE ===============================================================
%
%   figure; 
%
%   nX = [0:.001:0.05]; 
%   nF = forceten(nX);
%   
%   plot(nX, nF);
%
%   title ('dimensionless tendon force-length relationship');
%   xlabel('tendon length, n.u.');
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

% According to [1,2], the tendon model can be described as follows.
fl_ten = @(x) 0.10377*(exp(91*x)-1).*(x<=0.01516) + ...
    (37.526*x-0.26029).*(x>0.01516);

% Compute output tendon force. 
nF = fl_ten(nX); 

end % function