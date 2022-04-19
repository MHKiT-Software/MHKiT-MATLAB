function [heading, pitch, roll] = orient2euler(omat)
%     Calculate DOLfYN-defined euler angles from the orientation matrix.
% 
%     Parameters
%     ----------
%     omat : array
%       The orientation matrix
% 
%     Returns
%     -------
%     heading : |ndarray|
%       The heading angle. Heading is defined as the direction the x-axis 
%       points, positive clockwise from North (this is *opposite* the 
%       right-hand-rule around the Z-axis), range 0-360 degrees.
%     pitch : ndarray
%       The pitch angle (degrees). Pitch is positive when the x-axis 
%       pitches up (this is *opposite* the right-hand-rule around the 
%       Y-axis).
%     roll : ndarray
%       The roll angle (degrees). Roll is positive according to the 
%       right-hand-rule around the instrument's x-axis.

    if isstruct(omat)
        if isfield(omat,'orientmat')
            omat = omat.orientmat.data;
        end
    end
    
    heading = rad2deg(atan2(omat(:,:,1,1),omat(:,:,2,1)));
    heading = rem(heading,360);
    
    pitch = rad2deg(asin(omat(:,:,3,1)));
    
    roll = rad2deg(atan2(omat(:,:,3,2),omat(:,:,3,3)));
    
end

