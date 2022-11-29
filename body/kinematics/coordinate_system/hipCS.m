function [i_L, j_L, k_L, i_R, j_R, k_R] = hipCS(fle, fme, uPh_L, uPh_R, uPk_L, uPk_R)
 %HIPCS Set up hip coordinate system.
    %   
    %   [i_L, j_L, k_L, i_R, j_R, k_R] = hipCS(fle, fme, uPh_L, uPh_R, uPk_L, uPk_R)
    %
    %   INPUT ========================================================
    %
    %   fle (structure)
    %   Contains triaxial positions of right and left FLE markers in 
    %       global coordinate system. 
    %   Example: fle.L.x = [1], fle.L.y = [1], fle.L.z = [1]
    %            fle.R.x = [1], fle.R.y = [1], fle.R.z = [1]  
    %
    %   fme (structure)
    %   Contains triaxial positions of right and left FME markers in 
    %       global coordinate system. 
    %   Example: fme.L.x = [1], fme.L.y = [1], fme.L.z = [1]
    %            fme.R.x = [1], fme.R.y = [1], fme.R.z = [1]
    %
    %   uPh_L (numeric array)
    %   Left hip joint center
    %   Example: [1;1;1]
    %
    %   uPh_R (numeric array)
    %   Right hip joint center
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
    %   Unit vector poiting forward. 
    %   
    %   j_* (numeric array)
    %   Unit vector poiting in the medial-lateral direction. 
    %   
    %   k_* (numeric array)
    %   Unit vector poiting along the tibia. 
    %
    %   * corresponds to the left and right legs
    %
    %   AUTHOR =========================================================
    %
    %   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
    %
    %   ================================================================

    % Find unit vector pointing along the femur
    k_L = uPk_L - uPh_L; k_L = k_L/vecnorm(k_L); 
    k_R = uPk_R - uPh_R; k_R = k_R/vecnorm(k_R);
    
    % Define vector laying in the plane with k_*
    uV_L = [fle.L.x - fme.L.x; fle.L.y - fme.L.y; fle.L.z - fme.L.z];
    uV_R = [fle.R.x - fme.R.x; fle.R.y - fme.R.y; fle.R.z - fme.R.z];
    
    % Find unit vector pointing forward (check the direction of rotation)
    i_L = cross(k_L, uV_L); i_L = i_L/vecnorm(i_L);
    i_R = cross(uV_R, k_R); i_R = i_R/vecnorm(i_R);
    
    % Find unit vector pointing in medial-lateral direction
    j_L = cross(k_L, i_L); 
    j_R = cross(k_R, i_R); 

end

