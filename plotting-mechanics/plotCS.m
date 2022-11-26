function [] = plotCS(O,X,Y,Z)
%PLOTCS Plots coordinate system's axes and its origin. 
%
%   plotCS(O,X,Y,Z)
%   
%   INPUT ===========================================================  
%
%   X (numeric array)
%   Contains spatial coordinates of the vector X. 
%   Example: [1,1,1]
%
%   Y (numeric array)
%   Contains spatial coordinates of the vector Y. 
%   Example: [1,1,1]
%
%   Z (numeric array)
%   Contains spatial coordinates of the vector Z. 
%   Example: [1,1,1]
%
%   O (numeric array)
%   Contains spatial coordinates of the origin of coordinate system. 
%   Example: [1,1,1]
%
%   AUTHOR =========================================================
%
%   S.Bahdasariants, NEL, WVU
%
%   ================================================================

    % PLOT X. Plot vector X originating in O
    plot3([O(1),X(1)+O(1)], [O(2),X(2)+O(2)], [O(3),X(3)+O(3)], 'r',...
        'LineWidth', 2); hold on; 
    
    % PLOT Y. Plot vector Y originating in O
    plot3([O(1),Y(1)+O(1)], [O(2),Y(2)+O(2)], [O(3),Y(3)+O(3)], 'g',...
        'LineWidth', 2); hold on; 
    
    % PLOT Z. Plot vector Z originating in O
    plot3([O(1),Z(1)+O(1)], [O(2),Z(2)+O(2)], [O(3),Z(3)+O(3)], 'b',...
        'LineWidth', 2); hold on; 
    
    % PLOT O. Plot origin
    plot3([O(1), O(1)], [O(2), O(2)], [O(3), O(3)], '.');
    
    % Visualize results in 3D
    view(3); 

end

