function nJ = paraxt(nI, nVec, nM)
%PARAXT Parallel axis theorem. 
%
%   INPUT =================================================================
%   
%   nI (numeric array)
%   3*3 moments of inertia matrix defined with respect to the old reference
%       frame, kg*m^2. 
%   Example: I = eye(3)
%
%   nVec (numeric array)
%   Coordinates of the vector pointing from the origin of the old reference
%       frame to the origin of the new reference frame.
%   Example: [1, 1, 1]
%
%   nM (double)
%   Mass of the object, kg. 
%   Example: 1
%
%   OUTPUT ================================================================
%
%   nJ (numeric array)
%   3*3 moments of inertia matrix defined with respect to the new reference
%       frame, kg*m^2. 
%   
%   EXAMPLE ===============================================================
%
%   I = eye(3); nCOM = [1,1,1]; nM = 1; J = paraxt(I, nCOM, nM); 
%   
%   AUTHOR ================================================================
%   
%   S.Bahdasariants, NEL, WVU, sb0220@mix.wvu.edu
%
%   See also MAIN, SOLVEDYNAMICS, SIMSWING, EXTRACTMETA, SAVESIM,...
%   GETSWING, RUNSIM, GETKIN, SCALEANTHRO, SETMECHANICS, FRUSTUMINERT,...
%   SETCIRCUM
%
%   =======================================================================

    % Preallocate moment of inertia matrix
    nJ = zeros(size(nI)); 
    
    % For all rows
    for i = 1:size(nI,1)
        % For all columns
        for j = 1:size(nI,2)
            % The inertia matrix defined with respect to the new origin of 
            % the reference frame:
            nJ(i,j) = nI(i,j) + nM*(vecnorm(nVec)*kroneckerDelta(i,j) - ...
                nVec(i)*nVec(j)); 
        end % j
    end % i
end % paraxt

function bDelta = kroneckerDelta(i,j)
% KRONECKERDELTA Kronecker delta is a function of two variables, usually 
% just non-negative integers. The function is 1 if the variables are equal,
% and 0 otherwise. 

    if i==j
        bDelta = 1; 
    else
        bDelta = 0; 
    end

end % kroneckerDelta

