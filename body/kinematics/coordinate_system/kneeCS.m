function [i_L, j_L, k_L, i_R, j_R, k_R] = kneeCS(uPa_L, uPa_R, uPk_L, uPk_R, MedioLat)
 %KNEECS Set up knee coordinate system.
    %   
    %   [i_L, j_L, k_L, i_R, j_R, k_R] = kneeCS(uPa_L, uPa_R, uPk_L, uPk_R, MedioLat)
    %
    %   INPUT ========================================================
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
    %   MedioLat (structure)
    %   Contains vectors pointing in medio-lateral direction for the right
    %       and left hips.
    %   Example: MedioLat.L = [1,1,1]; MedioLat.R = [1,1,1]
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

    % Find unit vector pointing along the tibia
    k_L = uPa_L - uPk_L; k_L = k_L/vecnorm(k_L); 
    k_R = uPa_R - uPk_R; k_R = k_R/vecnorm(k_R); 
    
    % Find unit vector pointing forward (check the direction of rotation)
    %   reusing the vector from the hip, poiting in medio-lateral direction
    i_L = cross(MedioLat.L, k_L); i_L = i_L/vecnorm(i_L); 
    i_R = cross(MedioLat.R, k_R); i_R = i_R/vecnorm(i_R); 
    
    % Reuse vector from the hip poiting in medial-lateral direction
    j_L = cross(k_L, i_L); 
    j_R = cross(k_R, i_R); 
    
end

