function nC = varcov(nX,nY)
%VARCOV Variance-covariance matrix. 
%   Computes variance-covariance (a.k.a. covariance) matrix for two inputs.
%   The diagonal terms correspond to the variances and the other entries to
%   covariances. 
%
%   nC = varcov(nX,nY)
%
%   INPUT =================================================================
%
%   nX (numeric array)
%   First input signal. 
%   Example: [0:.001:0.05]; 
%
%   nY (numeric array)
%   Second input signal. 
%   Example: [0:.01:0.5]; 
%
%   OUTPUT ================================================================
%
%   nF (numeric array)
%   Variance-covariance matrix.  
%
%   EXAMPLE ===============================================================
%
%   nX = [0:.001:0.05]; 
%   nY = [0:.01:0.5]; 
%   nC = varcov(nX,nY)
% 
%   REFERENCES ============================================================
%
%   1. https://datascienceplus.com/understanding-the-covariance-matrix/
%
%   AUTHOR ================================================================
%   
%   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
%
%   =======================================================================

% Compute the number of samples in each signal. 
nSmpl_x = length(nX); nSmpl_y = length(nY); 

% If the number of samples is different, throw an error. 
if nSmpl_x ~= nSmpl_y
    error('The number of elements in x and y must match.'); 
end

% % Compute each entry of the variance-covariance matrix separately. 
% nC(1,1) = 1/(nSmpl_x-1)*(nX - mean(nX))'*(nX-mean(nX)); 
% nC(1,2) = 1/(nSmpl_x-1)*(nX - mean(nX))'*(nY-mean(nY)); 
% nC(2,1) = 1/(nSmpl_x-1)*(nY - mean(nY))'*(nX-mean(nX)); 
% nC(2,2) = 1/(nSmpl_x-1)*(nY - mean(nY))'*(nY-mean(nY)); 

nC = [nX(:)' - mean(nX); nY(:)' - mean(nY)]*...
     [nX(:)' - mean(nX); nY(:)' - mean(nY)]'/(nSmpl_x-1); 

end