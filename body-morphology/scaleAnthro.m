function [nLength, nDim1, nDim2, nMass, nCOM, nInertiaMom,...
    nIntertiaProd, nVec] = scaleAnthro(nHeight, nWeight, nAge, sSex)
    %SCALEANTHRO Creates a structure containing subject's anthropometry. 
    %
    %   INPUT==============================================================
    %
    %   nHeight (double)
    %   Height of the subject (m).
    %   Example: 1.74
    %
    %   nWeigth (double)
    %   Weigth of the subject (kg). 
    %   Example: 70
    %
    %   nAge (double)
    %   Age of the subject.
    %   Example: 50
    %
    %   sSex (string)
    %   Sex of the subject. 
    %   Example: "M"
    %   
    %   OUTPUT=============================================================
    %   
    %   nLength (structure)
    %   Contains the lengths of the segments.
    %
    %   nDim1 (structure)
    %   Contains the widths and radii of the segments. 
    %   
    %   nDim2 (structure)
    %   Contains the heights of the segments. 
    %
    %   nMass (structure)
    %   Contains masses of the segments. 
    %
    %   nCOM (structrue)
    %   Coontains coordinates of the center of mass for each segment. 
    %
    %   nInertiaMom (structure)
    %   Contains moments fo inertia for each segment. 
    %
    %   nIntertiaProd (structure)
    %   Contains inertia products for each segment. 
    %
    %   nVec (structure)
    %   Coordinates of the vector pointing from the origin of the 
    %   reference frame to the hinge point.
    %
    %   AUTHOR ============================================================
    %   
    %   S.Bahdasariants, NEL, WVU, sb0220@mix.wvu.edu
    %
    %   See also MAIN, SOLVEDYNAMICS, SIMSWING, EXTRACTMETA, SAVESIM,...
    %   GETSWING, RUNSIM, SETCIRCUM, GETKIN, SETMECHANICS, PARAXT,...
    %   FRUSTUMINERT
    %
    %   LITERATURE ========================================================
    %    
    %   1. Winter, David A. Biomechanics and motor control of human 
    %   movement. John Wiley & Sons, 2009, p. 83-86.
    % 
    %   2. Fryar, Cheryl D., Margaret D. Carroll, Qiuping Gu, Joseph Afful,
    %   and Cynthia L. Ogden. 2021. “Anthropometric Reference Data for 
    %   Children and Adults: United States, 2015-2018.”
    %
    %   3. Hanavan Jr, Ernest P. ‘A Mathematical Model of the Human Body’. 
    %   Air Force Aerospace Medical Research Lab Wright-patterson AFB OH,
    %   1964.
    %
    %   4. "Stature Estimation Based on Dimensions of the Bony Pelvis and
    %   Proximal Femur* - Giroux - 2008 - Journal of Forensic Sciences - 
    %   Wiley Online Library.” Accessed December 1, 2021.
    %
    %   5. The mathematical formulations, notations, and derivations can be 
    %   accessed through the following link: https://drive.google.com/file
    %   /d/1HVyCyX7EoZSMdfH6forwRxNf6-CP4Mux/view?usp=sharing
    %
    %   ===================================================================

    
    %% FOOT. Define dimensions and inertial properties of the foot. 
    
    % Define length of the foot [1]
    nLength.foot = 0.152*nHeight;

    % Define height of the foot [1] (Sphyrion Height)
    nSphyHeight = 0.039*nHeight;

    % Define bigger radius of the frustum of the right circular cone [3]
    nDim1.foot = 0.5*nSphyHeight; 

    % Find the smaller radius from the equation for the centroid of the
    % frustrum [3]:
   
    fun_foot = @(RR) 0.25*nLength.foot*(nDim1.foot^2 + 2*nDim1.foot*RR +...
        3*RR^2)/(nDim1.foot^2 + nDim1.foot*RR + RR^2) - 0.429*nLength.foot; 

    % Set options: do not print out fsolve result
    options = optimset('Display','off');

    nDim2.foot = fsolve(fun_foot,0.75*nDim1.foot, options);
    
    % Define mass of the foot, kg [1]
    nMass.foot = 0.0145*nWeight;

    % Coordinates of the centroid relative to the reference frame. The
    % center of the reference frame is set in the middle of the segment by
    % default [3]
    nCOM.foot = [0, 0, 0.5*nLength.foot - 0.429*nLength.foot];
    
    % Calculate tensor of moments of intertia with respect to the base of 
    % the object, kg*m^2:
    nI = frustumInert(nDim1.foot, nDim2.foot, nLength.foot, nMass.foot); 

    % Define vector pointing from the origin of the old reference frame 
    % (from the center of mass) to the hinge point [3,5]
    nNumerator = 0.25*nLength.foot*(nDim1.foot.^2 + 2*nDim1.foot*...
                 nDim2.foot + 3*nDim2.foot^2); 
    nDenominator = nDim1.foot.^2 + nDim1.foot*nDim2.foot + nDim2.foot.^2;
    nL3 = nNumerator/nDenominator; 
    nH1 = nLength.foot*nDim1.foot/(nDim1.foot-nDim2.foot);
    nL4 = nDim1.foot*(nH1 - nL3)/nH1; 
    nVec.foot = [-nL4, 0, nL3]; 

    % Use parrallel axis theorem to calculate the inertia with respect to
    % the hinge point:
    nJ = paraxt(nI, nVec.foot, nMass.foot); 

    % Moments of inertia
    nInertiaMom.foot = diag(nJ)';

    % Products of inertia
    nIntertiaProd.foot = [nJ(2,3), nJ(3,1), nJ(1,2)]; 
    
    %% THIGH. Define dimensions and inertial properties of the thighs. 
    
    % Define length of the thigh [1]
    nLength.thigh = 0.245*nHeight;
    
    % Fetch thigh circumference from age- and sex-based anthropometric
    % tables [2]:
    nCircum = setCircum(nAge,sSex, "thigh", "metaCircumference");
    
    % Define the bigger radius of the frustum of the right circular cone
    nDim1.thigh = nCircum/(2*pi);
    
    % Define the smaller radius of the frustrum of the right circ. cone
    fun_thigh = @(RR) 0.25*nLength.thigh*(nDim1.thigh^2 +...
        2*nDim1.thigh*RR + 3*RR^2)/(nDim1.thigh^2 + nDim1.thigh*RR +...
        RR^2) - 0.437*nLength.thigh;
    
    nDim2.thigh = fsolve(fun_thigh, 0.75*nDim1.thigh, options);
    
    % Define mass of the thigh
    nMass.thigh = 0.1*nWeight;
    
    % Coordinates of the centroid relative to the reference frame. The
    % center of the reference frame is set in the middle of the segment by
    % default
    nCOM.thigh = [0, 0, 0.437*nLength.thigh - 0.5*nLength.thigh];

    % Calculate tensor of moments of intertia with respect to the base of 
    % the object, kg*m^2:
    nI = frustumInert(nDim1.thigh, nDim2.thigh, nLength.thigh, nMass.thigh); 
    
    % Define vector pointing from the origin of the old reference frame 
    % (from the center of mass) to the hinge point [1, 3]
    nVec.thigh = [0, 0, -0.01*nHeight]; 

    % Use parrallel axis theorem to calculate the inertia with respect to
    % the hinge point:
    nJ = paraxt(nI, nVec.thigh, nMass.thigh); 

    % Moments of inertia
    nInertiaMom.thigh = diag(nJ)';

    % Products of inertia
    nIntertiaProd.thigh = [nJ(2,3), nJ(3,1), nJ(1,2)]; 
    
    %% SHANK. Define dimensions and inertial properties of the shanks.
    
    % Fetch shank circumference from age- and sex-based anthropometric
    % tables [2]:
    nLength.shank = 0.246*nHeight;
    
    % Shank circumference [2]
    nCircum = setCircum(nAge,sSex,"calf","metaCircumference");
    
    % Define the bigger radius of the frustum of the right circular cone[3]
    nDim1.shank = nCircum/(2*pi);
    
    % Define the smaller radius of the frustrum of the right circular cone
    fun_shank = @(RR) 0.25*nLength.shank*(nDim1.shank^2 +...
        2*nDim1.shank*RR + 3*RR^2)/(nDim1.shank^2 + nDim1.shank*RR +...
        RR^2) - 0.416*nLength.shank;
    
    nDim2.shank = fsolve(fun_shank, 0.75*nDim1.shank, options);
    
    % Define mass of the shank, kg
    nMass.shank = 0.0465*nWeight;
    
    % Coordinates of the centroid relative to the reference frame. The
    % center of the reference frame is set in the middle of the segment by
    % default
    nCOM.shank = [0, 0, 0.416*nLength.shank - 0.5*nLength.shank];
    
    % Calculate tensor of moments of intertia with respect to the base of 
    % the object, kg*m^2:
    nI = frustumInert(nDim1.shank, nDim2.shank, nLength.shank, nMass.shank); 
    
    % Define vector pointing from the origin of the old reference frame 
    % (from the center of mass) to the hinge point [1, 3]
    nVec.shank = [0, 0, 0]; 

    % Use parrallel axis theorem to calculate the inertia with respect to
    % the hinge point:
    nJ = paraxt(nI, nVec.shank, nMass.shank); 

    % Moments of inertia
    nInertiaMom.shank = diag(nJ)';

    % Products of inertia
    nIntertiaProd.shank = [nJ(2,3), nJ(3,1), nJ(1,2)]; 
    
    %% PELVIS. Define dimensions and inertial properties of the pelvis.       
    
    % Estimate length of the trunk from its mean circumference measured as 
    %   a function of age and sex in [2]
    nLength.trunk = setCircum(nAge,sSex,"waist","metaCircumference")/(2*pi);
    
    % Define trunk width [1]
    nDim1.trunk   = 0.191*nHeight;
    
    % The height of the trunk is estimated from the regression 
    %   equations[4]
    if sSex == "M"
        % Compute trunk height from the reg. eq. 
        nDim2.trunk = (nHeight - 1.04473)/3.262; 
    elseif sSex == "F"
        % Compute trunk height from the reg. eq. 
        nDim2.trunk = (nHeight - 0.92784)/3.480; 
    end % if
    
    % Mass of the pelvis is initialized with aggregated mass of hands,
    % forearms, upperarms, head, neck, and torso [1]
    nArm      = 0.050*nWeight;
    nHand     = 0.006*nWeight;
    nTrunk    = 0.497*nWeight;
    nHeadNeck = 0.081*nWeight;
    
    % Trunk mass (aggregated mass, see explanation above)
    nMass.trunk  = nTrunk + nHeadNeck + nArm*2 + nHand*2; 
    
    % Center of mass coordinates of the center of mass relative to the 
    %   block reference frame.
    nCOM.trunk = [0,0,0];
    
    % To place the right hip in the pelvis correctly, I introduce the 
    %   values of vertical and horizontal insertion of the femur into the 
    %   pelvis.
    nLength.pelvis = 0;            % vertical
    nDim1.pelvis   = nDim1.thigh;  % horizontal
    
end % function

