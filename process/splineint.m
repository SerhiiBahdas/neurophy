function Y = splineint(X,tOld, tNew)
    %SPLINEINT Resamples signal X on time tNew.
    %
    %   Y = splineint(X,tOld, tNew)
    %
    %   INPUT =============================================================
    %   
    %   x (numeric array)
    %   Input signal to be resampled. 
    %   Example: [1,2,3; 4,5,6]
    %
    %   tOld (numeric array)
    %   Time vector on which signal is defined. 
    %   Example: [0, 1, 2]
    %
    %   tNew (numeric array)
    %   Time vector on which signal needs to be interpolated. 
    %   Example: [0, 0.5, 1, 1.5, 2]
    %
    %   OUTPUT ============================================================
    %   
    %   X (numeric array)
    %   Resampled signal.
    %
    %   AUTHOR ============================================================
    %
    %   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
    %
    %   ===================================================================
    
    
    % Determine the size of the input matrix (rows and columns). The 
    %   greater dimention should contain time samples. The smaller 
    %    dimention is equal to the number of signals in the matrix. 
    
    [row, col] = size(X); 
    
    % Determine signal and sample dimentions 
    if row < col, nSig = row; else X = X'; nSig = col; end
    
    for iSignal = 1:nSig % loop through signals 
        
        % Assign resampled signal to the ouput matrix
        
        pp = spline(tOld, X(iSignal,:)); % piecewise polynomial structure
        Y(iSignal,:) = ppval(pp, tNew);  % resampled signal
    
    end
end

