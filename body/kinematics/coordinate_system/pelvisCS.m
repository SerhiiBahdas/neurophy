function [i, j, k] = pelvisCS(ias, ips)
    %PELVISCS Set up pelvis coordinate system.
    %   
    %   [i, j, k] = pelvisCS(ias, ips)
    %
    %   INPUT ========================================================
    %
    %   ias (structure)
    %   Contains triaxial positions of right and left IAS markers in 
    %       global coordinate system. 
    %   Example: ias.L.x = [1], ias.L.y = [1], ias.L.z = [1]
    %            ias.R.x = [1], ias.R.y = [1], ias.R.z = [1]
    %
    %   ips (structure)
    %   Contains triaxial positions of right and left IPS markers in 
    %       global coordinate system. 
    %   Example: ips.L.x = [1], ips.L.y = [1], ips.L.z = [1]
    %            ips.R.x = [1], ips.R.y = [1], ips.R.z = [1]
    %
    %   OUTPUT =======================================================
    %
    %   i (numeric array)
    %   Unit vector poiting in the direction of locomotion. 
    %   
    %   j (numeric array)
    %   Unit vector poiting in the medial-lateral direction. 
    %   
    %   k (numeric array)
    %   Unit vector poiting upward. 
    %
    %   AUTHOR =========================================================
    %
    %   S.Bahdasariants, NEL, WVU
    %
    %   ================================================================
    
    % Find a midpoint between right and left IAS markers
    iasm = [ias.L.x + ias.R.x; ias.L.y + ias.R.y; ias.L.z + ias.R.z]./2;
    
    % Find midpoint between right and left IPS markers
    ipsm = [ips.L.x + ips.R.x; ips.L.y + ips.R.y; ips.L.z + ips.R.z]./2;
      
    % Find a vector laying in the same plane with IPS and IAS
    v1 = iasm - ipsm; 
    
    % Find a vector poiting from the midpoint between IAS markers to the 
    %   left IAS marker
    v2 = [ias.R.x; ias.R.y; ias.R.z] - iasm; 
    
    % Find a unit vector pointing in medial-lateral direction
    j = v2/vecnorm(v2); 
    
    % Find a vector poiting down
    v3 = cross(v1, j); 
    
    % Find unit vector pointing down
    k = v3/vecnorm(v3); 
    
    % Find vector poiting in the direction of the movement. 
    i = cross(j,k); 
    
end

