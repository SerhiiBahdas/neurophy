function nF = forcelen(nX, varargin)
%FORCELEN Force-length relationship describes the dependence of the
% maximum isometric force on the length of the muscle [fiber, sarcomere]. 
%
%   nF = forcelen(nX)
%   nF = forcelen(nX, 'sMethod', 'Yakovenko')
%
%   INPUT =================================================================
%
%   nX (numeric array)
%   Length of the muscle, m. 
%   Example: [0:0.05:1]; 
%
%   [OPTIONAL INPUT]
%
%   sMethod (string)
%     
% 
% 
% The original Hill's approximation of the force-length relationship pro-
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


end