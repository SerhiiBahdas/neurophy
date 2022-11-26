function [] = vecplot(basis, uPo, uP1, uP2, bCSplot)
    %VECPLOT Plots a vector between two points in the input basis 
    %   coordinates.
    %
    %   vecplot(basis, uPo, uP1, uP2, bCSplot)
    %
    %   INPUT =========================================================
    %
    %   basis (structure)
    %   Contains orthonormal basis' vectors.
    %   Example: basis.ex = [1;0;0]; 
    %            basis.ey = [0;1;0]; 
    %            basis.ez = [0;0;1];
    %
    %   uPo (numeric array)
    %   Origin of the input basis in global coordinate system. 
    %   Example: [0;0;0]
    %
    %   uP1 (numeric array)
    %   Triaxial coordinates of the vector's start point in global 
    %       coordinate system.
    %   Example: [1;1;1]
    %
    %   uP2 (numeric array)
    %   Triaxial coordinates of the vector's end point in global 
    %       coordinate system.
    %   Example: [2;2;2] 
    %
    %   bCSplot (boolean)
    %   Plot coordinate system defined by input orthonormal basis. 
    %   Example: 1, 0
    %
    %   CODE EXAMPLE =================================================
    %
    %   % Create input basis structure 
    %   basis.ex = rotx(30)*[1;0;0];
    %   basis.ey = rotx(30)*[0;1;0];
    %   basis.ez = rotx(30)*[0;0;1];
    %   
    %   % Specify the origin of the basis (in global coordinate system)
    %   uPo = [2,2,2]';
    %   
    %   % Initial and terminal points of the vector (in global coordinate 
    %   system). Vector points from point 1 to point 2.
    %
    %   uP1 = [2,4,1]'; uP2 = [3,3,3]';
    %
    %   % Use VECPLOT to visualize the basis and the vector in basis' 
    %       coordinate system. 
    %
    %   vecplot(basis, uPo, uP1, uP2, 1)
    %
    %
    %   REFERENCES ===================================================
    %
    %   1. Valero-Cuevas, Francisco J. Fundamentals of neuromechanics. 
    %       Vol. 8. Berlin: Springer, 2016.
    %
    %   AUTHORSHIP ===================================================
    %
    %   Serhii Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
    %
    %   ==============================================================
    
    % DATA PREPROCESSING. Compute vector coordinates.
    
    % Calculate vector coordinates in global coordinate system 
    %   from two points
    uP = uP2 - uP1; 

    % Define global coordinate system vector basis
    u.ex = [1;0;0]'; u.ey = [0;1;0]'; u.ez = [0;0;1]';
    
    % Compute rotation matrix from global coordinate system to the CS 
    %   defined with input orthonormal vector basis (see Linear Algebra 
    %   textbook for kindergarteners, à la Sergiy, or see [1])
    
    bRu = [dot(basis.ex, u.ex), dot(basis.ex, u.ey), dot(basis.ex, u.ez); 
           dot(basis.ey, u.ex), dot(basis.ey, u.ey), dot(basis.ey, u.ez); 
           dot(basis.ez, u.ex), dot(basis.ez, u.ey), dot(basis.ez, u.ez)];
    
    % Compute vector coordinates in the input orthonormal basis
    bP = bRu*(uP - uPo); 

    % Plot the vector originating in O. 
    plot3([uPo(1),bP(1) + uPo(1)], [uPo(2),bP(2)+ uPo(1)],...
        [uPo(3),bP(3)+uPo(1)], 'k', 'LineWidth', 2);
    hold on;
    
    % PLOT O. Plot the origin. 
    plot3([uPo(1), uPo(1)], [uPo(2), uPo(2)], [uPo(3), uPo(3)], '.');
    hold on;
 
    if bCSplot == 1 % Plot coordinate system defined by input basis
        
        % PLOT X. Plot vector X originating in O. 
        plot3([uPo(1),basis.ex(1) + uPo(1)], [uPo(2),basis.ex(2) + uPo(2)],...
            [uPo(3),basis.ex(3) + uPo(3)], 'r', 'LineWidth', 2); 
        hold on; 

        % PLOT Y. Plot vector Y originating in O. 
        plot3([uPo(1),basis.ey(1) + uPo(1)], [uPo(2),basis.ey(2) + uPo(2)],...
            [uPo(3),basis.ey(3) + uPo(3)], 'g', 'LineWidth', 2); 
        hold on; 

        % PLOT Z. Plot vector Z originating in O. 
        plot3([uPo(1),basis.ez(1) + uPo(1)], [uPo(2),basis.ez(2) + uPo(2)],...
            [uPo(3),basis.ez(3) + uPo(3)], 'b', 'LineWidth', 2); 
        hold on; 

        % Visualize results in 3D. 
        view(3); axis equal; grid on; 
        
    end % end if bCSplot
end % function