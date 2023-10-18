function y = singleTwitch_SF(x, bPlot)
% singleTwitch_SF Computes the y values based on the Gaussian Curve Fit.
%
%   y = singleTwitch_SF(x)
%
%   INPUT =================================================================
%
%   x (nSamples x 1 numeric array)
%   X-axis values corresponding to the dataset.
%   Example: x = linspace(0, 10, 1000);
%
%   bPlot (boolean)
%   Flag for plotting. 
%   Example: 1
%
%   OUTPUT ================================================================
%
%   y (nSamples x 1 numeric array)
%   Resultant y values based on the Gaussian Curve Fit.
%
%   The Gaussian Curve Fit is defined as:
%   f(x) = a1*exp(-((x-b1)/c1)^2) + a2*exp(-((x-b2)/c2)^2) + a3*exp(-((x-b3)/c3)^2) 
%        + a4*exp(-((x-b4)/c4)^2) + a5*exp(-((x-b5)/c5)^2)
%
%   Coefficients for Gaussian Curve Fit are as follows:
%   Coefficients (with 95% confidence bounds):
%       a1 = 1 (default value, as not provided)
%       b1 = -82.6425
%       c1 = 62.1849
%       a2 = -0.1072
%       b2 = 89.5301
%       c2 = 15.5481
%       a3 = -3.1335
%       b3 = 17.1072
%       c3 = 43.0769
%       a4 = 2.8766e+04
%       b4 = -585.3943
%       c4 = 215.5835
%       a5 = -0.7172
%       b5 = 12.0961
%       c5 = 14.8576
%
%   AUTHOR ================================================================
%
%   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
%
%   =========================================================================

% Coefficients
a1 = -87.7936; 
b1 = -82.6425;
c1 = 62.1849;
a2 = -0.1072;
b2 = 89.5301;
c2 = 15.5481;
a3 = -3.1335;
b3 = 17.1072;
c3 = 43.0769;
a4 = 2.8766e+04;
b4 = -585.3943;
c4 = 215.5835;
a5 = -0.7172;
b5 = 12.0961;
c5 = 14.8576;

% Convert time to ms. 
x = x*1000; 

% Gravity Constant
g = 9.81; 

y = (a1.*exp(-((x-b1)./c1).^2) + a2.*exp(-((x-b2)./c2).^2) + a3.*exp(-((x-b3)./c3).^2) ...
    + a4.*exp(-((x-b4)./c4).^2) + a5.*exp(-((x-b5)./c5).^2))/g/1000;

% Convert time to ms. 
x = x/1000; 

% Plotting
if bPlot 

    figure; 
    plot(x,y, 'LineWidth', 1)
    xlabel('Time, s')
    ylabel('Amplitude (slow fatiguable), N')
    title('Change in muscle tension in response to a single action potential')

end

end
