function out = calc_beam_orientmat(theta, convex, degrees)
%%%%%%%%%%%%%%%%%%%%
%    Calculate the rotation matrix from beam coordinates to
%    instrument head coordinates for an RDI ADCP.
%     
% Parameters
% ------------
%    theta : is the angle of the heads (usually 20 or 30 degrees)
%
%    convex : is a flag for convex or concave head configuration.
%
%    degrees : is a flag which specifies whether theta is in degrees
%        or radians (default: degrees=True)
%
% Returns
% ---------
%     out: array 
%         orientation matrix
%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if degrees
        theta = deg2rad(theta);
    end

    if convex == 0 || convex == -1
        c = -1;
    else
        c = 1;
    end

    a = 1 / (2. * sin(theta));
    b = 1 / (4. * cos(theta));
    d = a / (2.^0.5);

    out = [c*a, -c*a,    0,   0;
             0,    0, -c*a, c*a;
             b,    b,    b,   b;
             d,    d,    d,  -d];
    
end

