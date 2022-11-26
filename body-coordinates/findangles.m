function [yRot, xRot, zRot] = findangles(proximal, distal)
    %FINDANGLES Finds joint angles from the rotation between coordinate
    %   systems. 
    %
    %   [yRot, xRot, zRot] = findangles(proximal, distal)
    %
    %   INPUT =========================================================
    %
    %   proximal (structure)
    %   Orthonormal basis of the proximal segment/joint. 
    %   Example: proximal.i = [1;0;0]
    %            proximal.j = [0;1;0]
    %            proximal.k = [0;0;1]
    %
    %   distal (structure)
    %   Orthonormal basis of the distal segment/joint. 
    %   Example: distal.i = [1;0;0]
    %            distal.j = [0;1;0]
    %            distal.k = [0;0;1]
    %
    %   OUTPUT =======================================================
    %
    %   xRot (double)
    %   Rotation (an angle) around x axis, deg. 
    %
    %   yRot (double)
    %   Rotation (an angle) around y axis, deg. 
    %
    %   zRot (double)
    %   Rotation (an angle) around z axis, deg. 
    %
    %   AUTHOR =========================================================
    %
    %   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
    %
    %   ================================================================
    
    % Find rotation from proximal to distal (see "Fundamentals of 
    % Neuromechanics", Francisco J.Valero-Cuevas, p.184)
    dRp = [dot(distal.i, proximal.i), dot(distal.i, proximal.j), dot(distal.i, proximal.k);...
           dot(distal.j, proximal.i), dot(distal.j, proximal.j), dot(distal.j, proximal.k); 
           dot(distal.k, proximal.i), dot(distal.k, proximal.j), dot(distal.k, proximal.k)];

    % Find rotation around each axis
    nAng = rotm2eul(dRp, "ZYX"); 
    
    % X ROTATION. Find rotation x axis (Abduction-Adduction)
    xRot = rad2deg(nAng(3)); 
    
    % Y ROTATION. Find rotation around y axis (Flexion-Extension)
    yRot = rad2deg(nAng(2)); 
    
    % Z ROTATION. Find rotation around z axis (Int-Ext Rotation)
    zRot = rad2deg(nAng(1));    
       
%     IF THE TOOLBOX USED ABOVE IS NOT AVAILABLE, USE: 
%
%     % Y ROTATION. Find rotation around y axis 
%     yRot = asin(pRd(3,1)); 
% 
%     % Z ROTATION. Find rotation around z axis 
%     zRot = asin(-pRd(2,1)/cos(yRot)); 
% 
%     % Determine if which quadrant the angle is 
%     if (-pRd(2,1)/cos(yRot)) > 0
% 
%         zRot = pi - zRot;
% 
%     elseif (-pRd(2,1)/cos(yRot)) > pi
% 
%         zRot = zRot - 2*pi; 
% 
%     end % end if 
% 
%     % X ROTATION. Find rotation x axis
%     xRot = asin(-pRd(3,2)/cos(yRot));
% 
%     % Determine if which quadrant the angle is 
%     if (-pRd(3,2)/cos(yRot)) > 0
% 
%         zRot = pi - zRot;
% 
%     elseif (-pRd(3,2)/cos(yRot)) > pi
% 
%         zRot = zRot - 2*pi; 
% 
%     end % end if 

end % function


