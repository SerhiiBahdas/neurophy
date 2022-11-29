% Specify variables. 
syms F K s

% Describe Ib firing rate as a function of s and F. 
fib = K*F*(s+0.15)*(s+1.5)*(s+16)/(s+0.2)/(s+2)/(s+37); 

% Apply inverse Laplas transformation to Describe Ib 
% firing rate as a function of t and F. 
eq = ilaplace(fib);

% Visualize equation. 
pretty(eq)