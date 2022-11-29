function [signal] = getderiv(nSignal, nT, nTNew)
    %GETDERIV resamples input signal and computes its two derivatives using
    %    cubic splite interpolation. 
    %
    %   [signal] = getderiv(nSignal, nT, nTNew)
    %
    %   INPUT ==========================================================
    %
    %   nSignal (numeric array)
    %   Input signal (positions, angles). 
    %   Example: [1,2,3]
    %
    %   nTime (numeric array)
    %   Old time vector.
    %   Example: [1,2,3]
    %
    %   nTimeNew (numeric array)
    %   New time vector.
    %   Example: [1, 1.5, 2, 2.5, 3]
    %
    %   OUTPUT =========================================================
    %
    %   signal (structure)
    %   Contains resampled positions, computed velocities and 
    %       acceelrations, and time.
    %
    %   AUTHOR =========================================================
    %
    %   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
    %
    %   ================================================================
    
    pp1  = spline(nT, nSignal); % joint angles coeffs
    pp2  = fnder(pp1,1);        % angular velocity coeffs
    pp3  = fnder(pp2,1);        % angular acceleration coeffs
    
    signal.pos = ppval(pp1,nTNew);
    signal.vel = ppval(pp2,nTNew);
    signal.acc = ppval(pp3,nTNew);
    signal.time   = nTNew; 

end

