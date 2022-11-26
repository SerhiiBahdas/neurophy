function bEqual = approxequal(nVal1, nVal2, nTol)
%APPROXEQUAL Finds if one value is approximately equal to another value 
% with specified tolerance.
%
%   bEqual = approxequal(nVal1, nVal2, nTol)
%
%   INPUT ==========================================================
%
%   nVal1 (double)
%   First value to be compared.
%   Example: 1
%
%   nVal2 (double)
%   Second value to be compared. 
%   Example: 1.05
%
%   nTol (double)
%   Tollerance. 
%   Example: 0.05
%   
%   OUTPUT =========================================================
%
%   bEqual (boolean)
%   1 if the values are equal, 0 if the values are not equal.
%
%   AUTHOR =========================================================
%
%   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
%
%   ================================================================

bEqual = (nVal1 <= (nVal2 + nTol)) & (nVal1 >= (nVal2 - nTol)); 

end