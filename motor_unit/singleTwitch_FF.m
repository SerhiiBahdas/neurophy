function y = singleTwitch_FF(x, bPlot)
%SINGLETWITCH_FF Gaussian Fit of Muscle Tension Response to Action Potential
%
%   y = FF(x)
%
%   This function performs a Gaussian fit of the change in muscle tension in response to a single action potential 
%   using a Gaussian Curve Fit (gauss8). The function returns the fitted y values for the given x values.
%
%   INPUT ==================================================================
%
%   x (numeric array)
%   The x-axis values where the Gaussian fit is to be evaluated.
%   Example: x = linspace(0, 10, 100);
%
%   bPlot (boolean)
%   Flag for plotting. 
%   Example: 1
%
%   OUTPUT =================================================================
%
%   y (numeric array)
%   The y-axis values of the Gaussian fit corresponding to the given x values.
%
%   COEFFICIENTS ===========================================================
%
%   The Gaussian fit is performed using predefined coefficients obtained from a previous curve fitting process. 
%   The coefficients and their 95% confidence bounds are as follows:
%
%   Coefficients:
%   b1 = 39.8878; c1 = 17.3750; a2 = 36.2756; b2 = 62.9144; c2 = 27.0419;
%   a3 = 12.0962; b3 = 28.8204; c3 = 10.1693; a4 = 12.7101; b4 = 93.9349;
%   c4 = 24.6594; a5 = 6.0025;  b5 = 22.2440; c5 = 5.1468;  a6 = 25.3883;
%   b6 = 114.7191; c6 = 40.6630; a7 = -0.1038; b7 = 121.3560; c7 = 0.8919;
%   a8 = 16.8146;  b8 = 162.2733; c8 = 54.0776;
%
%   95% Confidence Bounds:
%   [b1, c1, a2, b2, c2, a3, b3, c3, a4, b4, c4, a5, b5, c5, a6, b6, c6, a7, b7, c7, a8, b8, c8]
%   Lower bounds: [31.8640, 6.5499, 3.8170, 41.8713, 7.8434, -6.6136, 25.7772, 4.8984, -78.0588, 57.3962, -12.3649, 1.1765, 21.8722, 3.7643, -39.6886, 53.2356, -14.0563, -4.0378e+12, -1.6968e+13, -4.1589e+12, -2.0122, 116.0163, 35.8167]
%   Upper bounds: [47.9116, 28.2002, 68.7341, 83.9576, 46.2403, 30.8060, 31.8636, 15.4402, 103.4789, 130.4736, 61.6837, 10.8284, 22.6158, 6.5293, 90.4651, 176.2025, 95.3824, 4.0378e+12, 1.6968e+13, 4.1589e+12, 35.6415, 208.5304, 72.3384]
%
%   GOODNESS OF FIT ========================================================
%
%   R-square    0.9999 
%   DFE         45.0000
%   Adj R-sq    0.9998 
%   RMSE        0.2067 
%
%   EXAMPLE =================================================================
%
%   % Define x values
%   x = linspace(0, 10, 100);
%   % Call FF function
%   y = FF(x);
%   % Plot the Gaussian fit
%   figure;
%   plot(x, y);
%   title('Gaussian Fit of Muscle Tension Response to Action Potential');
%
%   AUTHOR ==================================================================
%
%   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
%
%   ========================================================================

% Found Coefficients:
a1 = 21.2924; 
b1 = 39.8878;
c1 = 17.3750;
a2 = 36.2756; 
b2 = 62.9144;
c2 = 27.0419; 
a3 = 12.0962; 
b3 = 28.8204; 
c3 = 10.1693;
a4 = 12.7101;
b4 = 93.9349;
c4 = 24.6594;
a5 = 6.0025;
b5 = 22.2440;
c5 = 5.1468;
a6 = 25.3883; 
b6 = 114.7191; 
c6 = 40.6630;
a7 = -0.1038;
b7 = 121.3560;
c7 = 0.8919;
a8 = 16.8146;
b8 = 162.2733;
c8 = 54.0776;

% Convert time to ms. 
x = x*1000; 

% Gravity Constant
g = 9.81; 

% Gaussian Fit Equation
y = (a1*exp(-((x-b1)./c1).^2) + a2*exp(-((x-b2)./c2).^2) + a3*exp(-((x-b3)./c3).^2) + a4*exp(-((x-b4)./c4).^2) + a5*exp(-((x-b5)./c5).^2) + a6*exp(-((x-b6)./c6).^2) + a7*exp(-((x-b7)./c7).^2) + a8*exp(-((x-b8)./c8).^2))/g/1000;

% Convert time back
x = x/1000; 

% Plotting
if bPlot 

    figure; 
    plot(x,y, 'LineWidth', 1)
    xlabel('Time, s')
    ylabel('Amplitude (fast fatiguable), N')
    title('Change in muscle tension in response to a single action potential')

end
end
