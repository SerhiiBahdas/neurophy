function y = singleTwitch_FFR(x, bPlot)
% singleTwitch_FFR Computes the y values based on the Gaussian Curve Fit.
%
%   y = singleTwitch_FFR(x)
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
%       b1 = 107.5836
%       c1 = 64.1433
%       a2 = 7.8410
%       b2 = 38.8494
%       c2 = 50.0925
%       a3 = -0.6838
%       b3 = 14.8590
%       c3 = 12.7950
%       a4 = -6.7547
%       b4 = -53.9542
%       c4 = 86.6767
%       a5 = 1.3210
%       b5 = 184.9326
%       c5 = 44.0938
%
%   AUTHOR ================================================================
%
%   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
%
%   =========================================================================

% Coefficients
a1 = 7.4817; 
b1 = 107.5836;
c1 = 64.1433;
a2 = 7.8410;
b2 = 38.8494;
c2 = 50.0925;
a3 = -0.6838;
b3 = 14.8590;
c3 = 12.7950;
a4 = -6.7547;
b4 = -53.9542;
c4 = 86.6767;
a5 = 1.3210;
b5 = 184.9326;
c5 = 44.0938;

% Convert time to ms. 
x = x*1000; 

% Gravity Constant
g = 9.81; 

% Computing y values based on the Gaussian Curve Fit
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
