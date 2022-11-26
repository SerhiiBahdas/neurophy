function nLL = leglen(ftc,fcc)
    %LEGLEN calculates the length of a leg.
    %
    %   INPUT ============================================================
    %
    %   ftc (structure)
    %   Contains triaxial positions of right and left FTC markers in 
    %       global coordinate system. 
    %   Example: ftc.L.x = [1], ftc.L.y = [1], ftc.L.z = [1]
    %            ftc.R.x = [1], ftc.R.y = [1], ftc.R.z = [1]
    %
    %   fcc (structure)
    %   Contains triaxial positions of right and left FCC markers in 
    %       global coordinate system. 
    %   Example: fcc.L.x = [1], fcc.L.y = [1], fcc.L.z = [1]
    %            fcc.R.x = [1], fcc.R.y = [1], fcc.R.z = [1]
    %
    %   OUTPUT ===========================================================
    %
    %   nLL (double)
    %   Leg length, m.
    %   Example: 0.8
    %
    %   AUTHOR =========================================================
    %
    %   S.Bahdasariants, NEL, WVU
    %
    %   ================================================================
    
    % Create function handle calculating euclidean distance (vector length)
    euclidist = @(a,b) sqrt(sum((a-b).^2));
    
    % Create vectors for greater tronchanter coordinates
    nFTC_L = [ftc.L.x, ftc.L.y, ftc.L.z];
    nFTC_R = [ftc.R.x, ftc.R.y, ftc.R.z];
    
    % Create vectors for posterious calcaneous coordinates
    nFCC_L = [fcc.L.x, fcc.L.y, fcc.L.z]; 
    nFCC_R = [fcc.R.x, fcc.R.y, fcc.R.z]; 
    
    % Calculate left and right leg lengths 
    nLL_L = euclidist(nFTC_L, nFCC_L); 
    nLL_R = euclidist(nFTC_R, nFCC_R); 
    
    % Mean leg length
    nLL = mean([nLL_L, nLL_R]); 
    
end

