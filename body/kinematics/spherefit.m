function [x0,y0,z0,r] = spherefit(X,Y,Z,varargin)
    %SPHEREFIT Finds center and radius of the sphere from the experimental
    %   data using least-squares method.
    %
    %   Inputs------------------------------------------------------------
    %   
    %   X (numeric array) 
    %   X-coordinate data recorded for the spherical movement of the 
    %   subject.
    %   
    %   Y (numeric array)
    %   Y-coordinate data recorded for the spherical movement of the 
    %   subject. 
    %
    %   Z (numeric array)
    %   Z-coordinate data recorded for the spherical movement of the 
    %   subject. 
    %
    %   tol (double)
    %   Optional. Specifies a tolerance for the least-squares method. 
    %   The default tolerance is 1e-6.
    %
    %   Outputs-----------------------------------------------------------
    %   
    %   x0 (double)
    %   x coordinate of the sphere's center 
    %   
    %   y0 (double)
    %   y coordinate of the sphere's center 
    %
    %   z0 (double
    %   z coordinate of the sphere's center 
    %   
    %   r (double)
    %   radius of the sphere
    %
    %   Method explanation------------------------------------------------ 
    %
    %   The equation of the sphere is: 
    %       (x-x0)^2 + (y-y0)^2 + (z-z0)^2 = r^2.
    %
    %   If exapnded, the equation can be written as follows:
    %       x^2 + y^2 + z^2 = 2*x*x0 + 2*y*y0 + 2*z*z0 + r^2 - x0^2 
    %           - y0^2 - z0^2
    %
    %   When written in the matrix form to include n experimental points, 
    %     the equation can be represented as a system of linear equations:
    %
    %   b = Ac, where A, b, and c are represented with  
    %
    %   b = [x(i)^2   + y(i)^2   + z(i)^2,
    %        x(i+1)^2 + y(i+1)^2 + z(i+1)^2,
    %        ...        ...        ...
    %        x(n)^2   + y(n)^2   + z(n)^2]
    %
    %   A = [ 2*x(i),   2*y(i),   2*z(i),   1;
    %         2*x(i+1), 2*y(i+1), 2*z(i+1), 1;
    %           ...         ...     ...     ...
    %        2*x(n),   2*y(n),   2*z(n),   1]
    % 
    %   c = [x0, y0, z0, r^2 - x0^2 - y0^2 - z0^2].
    %
    %   The task of determining the x0, y0, z0, and r can now be posed as
    %     an optimization problem of finding c. The optimization is done 
    %     with least-squares method. 
    %   
    %   See also LSQR
    
    
    % DEFINE VARIABLES. Define matrix A, and vector b, according to the 
    %   explanation given above. 
    
     b = X.*X + Y.*Y + Z.*Z;
     
     if ~iscolumn(b) % check if b is a column vector
         b = b';
     end
     
     A = zeros(length(X),4); % preallocate memory for the matrix
     
     A(:,1) = 2*X;  A(:,2) = 2*Y;   A(:,3) = 2*Z;   A(:,4) = 1;
     
     % OPTIMIZE. Solve the optimization problem with least-squares method. 
     %  If tollegrance value is provided, feed it into lsqr function. 
     
     switch nargin 
         case 3
            c = lsqr(A,b);
         case 4
            c = lsqr(A,b,varargin{1});
     end
     
     % OUTPUT. Caculate output variables. See METHOD EXPLANATION in 
     %   function's help section. 
     
     x0 = c(1); y0 = c(2);  z0 = c(3); % center of the rotation
     r = sqrt(c(4) + x0*x0 + y0*y0 + z0*z0); % radius of the sphere
     
end

