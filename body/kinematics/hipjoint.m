function hipCoord = hipjoint(nLL, ias, ftc)
    %HIPJOINT Finds hip joint centers in the pelvis coordinate system. 
    %
    %   hipCoord = hipjoint(nLL, ias, ftc)
    %
    %   INPUT =======================================================
    %   
    %   nLL (double)
    %   Leg length, m.
    %   Example: 0.8
    %
    %   ias (structure)
    %   Contains triaxial positions of right and left IAS markers in 
    %       global coordinate system. 
    %   Example: ias.L.x = [1], ias.L.y = [1], ias.L.z = [1]
    %            ias.R.x = [1], ias.R.y = [1], ias.R.z = [1]  
    %
    %   ftc (structure)
    %   Contains triaxial positions of right and left FTC markers in 
    %       global coordinate system. 
    %   Example: ftc.L.x = [1], ftc.L.y = [1], ftc.L.z = [1]
    %            ftc.R.x = [1], ftc.R.y = [1], ftc.R.z = [1]
    %
    %   OUTPUT =====================================================
    %
    %   hipCoord (structure)
    %   Hips joint centers in the pelvis coordinate system.
    %
    %   REFERENCES =================================================
    %
    %   1. Davis III, Roy B., Sylvia Ounpuu, Dennis Tyburski, and James R. 
    %   Gage. "A gait analysis data collection and reduction technique." 
    %   Human movement science 10, no. 5 (1991): 575-587.
    %
    %   2. Lencioni, Tiziana, Ilaria Carpinella, Marco Rabuffetti, Alberto 
    %   Marzegan, and Maurizio Ferrarin. "Human kinematic, kinetic and EMG
    %   data during different walking and stair ascending and descending 
    %   tasks." Scientific data 6, no. 1 (2019): 1-10.
    %
    %   AUTHOR =========================================================
    %
    %   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
    %
    %   ================================================================
    
    % Use constants from [1, p.583] to determine hip joint location
    nC = 0.115*nLL - 0.0153;  nT = deg2rad(28.4); nB = deg2rad(18);
    
    % Calculate distance between the landmarks in coronal plane
    nYdist = abs(ias.L.y - ias.R.y);
    
    % Specify the radius of the markers used (in meters) [2, p.2]
    nRmark = 6e-3; 
    
    % Estimate the distance (in meters) from the anterior superior iliac 
    %   spine to the great tronchanter in the saggital plane
    
    nXdist_R = abs(ias.R.x - ftc.R.x); 
    nXdist_L = abs(ias.L.x - ftc.L.x); 
    
    nXdist = mean([nXdist_L, nXdist_R]); 
    
    % The location (in meters) of the hip joint center in pelvis 
    %   coordinates relative to the origin of the pelvis embedded 
    %   coordinate system is defined as:
   
    Xh = (-nXdist - nRmark)*cos(nB) + nC*cos(nT)*sin(nB); 
    Yh = nC*sin(nT) - nYdist*0.5; 
    Zh = (-nYdist - nRmark)*sin(nB) - nC*cos(nT)*cos(nB); 
    
    % Hips joint centers in the pelvis coordinate system
    hipCoord.R = [Xh, -Yh, -Zh]; 
    hipCoord.L = [Xh,  Yh, -Zh]; 
   
end


