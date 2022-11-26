function [i_L, j_L, k_L, i_R, j_R, k_R] = ankleCS(fm5, fm1, uPa_L, uPa_R, uPk_L, uPk_R)
    %ANKLECS Set up ankle coordinate system.
    %   
    %   ankleCS(fm5, fm1, uPa_L, uPa_R, uPk_L, uPk_R)
    %
    %   INPUT ========================================================
    %
    %   fm5 (structure)
    %   Contains triaxial positions of right and left FM5 markers in 
    %       global coordinate system. 
    %   Example: fm5.L.x = [1], fm5.L.y = [1], fm5.L.z = [1]
    %            fm5.R.x = [1], fm5.R.y = [1], fm5.R.z = [1]
    %
    %   fm1 (structure)
    %   Contains triaxial positions of right and left FM1 markers in 
    %       global coordinate system. 
    %   Example: fm1.L.x = [1], fm1.L.y = [1], fm1.L.z = [1]
    %            fm1.R.x = [1], fm1.R.y = [1], fm1.R.z = [1]
    %
    %   uPa_L (numeric array)
    %   Left ankle joint center
    %   Example: [1;1;1]
    %
    %   uPa_R (numeric array)
    %   Right ankle joint center
    %   Example: [1;1;1]
    %
    %   uPk_L (numeric array)
    %   Left knee joint center
    %   Example: [1;1;1]
    %
    %   uPk_R (numeric array)
    %   Right knee joint center
    %   Example: [1;1;1]
    %
    %   OUTPUT =======================================================
    %
    %   i_* (numeric array)
    %   Unit vector poiting upward. 
    %   
    %   j_* (numeric array)
    %   Unit vector poiting in the medial-lateral direction. 
    %   
    %   k_* (numeric array)
    %   Unit vector poiting along the foot. 
    %
    %   * corresponds to the left and right legs
    %
    %   AUTHOR =========================================================
    %
    %   S.Bahdasariants, NEL, WVU
    %
    %   ================================================================
    
    % Find a midpoint between FM1 and FM5 markers
    uFM_L = [fm5.L.x + fm1.L.x; fm5.L.y + fm1.L.y; fm5.L.z + fm1.L.z]./2;
    uFM_R = [fm5.R.x + fm1.R.x; fm5.R.y + fm1.R.y; fm5.R.z + fm1.R.z]./2;
    
    % Find a unit vector pointing along the foot 
    k_L = uFM_L - uPa_L; k_L = k_L/vecnorm(k_L); 
    k_R = uFM_R - uPa_R; k_R = k_R/vecnorm(k_R); 
    
    % Find a vector pointing upward
    uV_L = uPk_L - uPa_L; 
    uV_R = uPk_R - uPa_R; 
    
    % Find unit vector pointing in medial-lateral direction
    j_L = cross(k_L, uV_L); j_L = j_L/vecnorm(j_L); 
    j_R = cross(k_R, uV_R); j_R = j_R/vecnorm(j_R); 
    
    % Find unit vector pointing upward
    i_L = cross(j_L, k_L); 
    i_R = cross(j_R, k_R); 
    
end

