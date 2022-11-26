function nI = frustumInert(nR1, nR2, nH, nM)
%FRUSTUMINERT Calculates moment of inertia matrix of the right circular
%   cone, defined with respect to the reference frame with the origin 
%   placed at the base of the frustum. 
%
%   nI = frustumInert(nR1, nR2, nH, nM)
%
%   The mathematical formulations, notations, and derivations can be 
%   accessed through the following link: https://drive.google.com/file/d/1
%   HVyCyX7EoZSMdfH6forwRxNf6-CP4Mux/view?usp=sharing
%
%   INPUT =================================================================
%
%   nR1 (double)
%   Bigger radius of the frustum, m. 
%   Example: 0.05
%
%   nR2 (double)
%   Smaller radius of the frustum, m.
%   Example: 0.03
%
%   nH (double)
%   Height of the frustum, m. 
%   Example: 0.1
%
%   nM (double)
%   Mass of the frustum, kg. 
%   Example: 1
%
%   OUTPUT ================================================================
%
%   nI (numeric array)
%   3*3 matrix of the moments of inertia. 
%
%   EXAMPLE ===============================================================
%
%   nR1 = 0.05; nR2 = 0.03; nH = 0.1; nM = 1; 
%   nI = frustumInert(nR1, nR2, nH, nM); 
% 
%   AUTHOR ================================================================
%   
%   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
%
%   =======================================================================

%% FIND ADDITIONAL PARAMETERS. Find heights and masses of the circular 
% cones forming the frustum. 

% Find heights of the two circular cones, m
nH1 = nH*nR1/(nR1-nR2); 
nH2 = nH*nR2/(nR1-nR2); 

% Define distance to the centers of mass for both circular cones, m
nL1 = nH1/4; 
nL2 = nH2/4; 

% Find volume of the bigger cone, m^3
nV1 = 1/3*pi*nR1*nR1*nH1; 
% Find volume of the smaller cone, m^3
nV2 = 1/3*pi*nR2*nR2*nH2; 
% Find volume of the frustum, m^3
nV = nV1 - nV2; 

% Under assumption of the equal density of the frustum along the three
% axes, the density can be calculated as (kg/m^3):
nRo = nM/nV; 

% Thus, the masses of the two cones are equal (kg):
nM1 = nRo*nV1; 
nM2 = nRo*nV2; 

% To find moment of intertia matrix for the frustum, we first find inertias
% of the two cones used to form frustum. All notations used can be accessed
% through the link provided in the description of the function.

% Bigger cone, kg*m^2
nIc   = 3/20*nM1*(nR1*nR1 + nH1*nH1/4); 
nIxib = nIc + nM1*nL1*nL1; 
nIzib = 0.3*nM1*nR1*nR1; 

% Smaller cone, kg*m^2
nIb   = 3/20*nM2*(nR2*nR2 + nH2*nH2/4); 
nIxis = nIb + nM2*(nL2 + nH)*(nL2 + nH); 
nIzis = 0.3*nM2*nR2*nR2; 

% Frustum, kg*m^2
nIx = nIxib - nIxis; 
nIy = nIx; 
nIz = nIzib - nIzis; 

% The matrix of moments of inertia, kg*m^2
nI = [nIx, 0, 0; 0, nIy, 0; 0, 0, nIz]; 

end % frustumInert









