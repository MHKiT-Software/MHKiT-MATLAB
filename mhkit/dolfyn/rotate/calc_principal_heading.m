function p_heading = calc_principal_heading(vel, tidal_mode)
% Compute the principal angle of the horizontal velocity.
%
%    Parameters
%    ----------
%    vel : array 
%      The 2D or 3D velocity array 
%    tidal_mode : bool 
%
%    Returns
%    -------
%    p_heading : float or ndarray
%      The principal heading in degrees clockwise from North.
%
%    Notes
%    -----
%    When tidal_mode=True, this tool calculates the heading that is
%    aligned with the bidirectional flow. It does so following these
%    steps:
%      1. rotates vectors with negative velocity by 180 degrees
%      2. then doubles those angles to make a complete circle again
%      3. computes a mean direction from this, and halves that angle
%         (to undo the doubled-angles in step 2)
%      4. The returned angle is forced to be between 0 and 180. So, you
%         may need to add 180 to this if you want your positive
%         direction to be in the western-half of the plane.
%
%    Otherwise, this function simply computes the average direction
%    using a vector method.

dt = vel(:,:,1) + vel(:,:,2) * 1i;

if tidal_mode
    % Flip all vectors that are below the x-axis
    dt(imag(dt) <= 0) = dt(imag(dt) <= 0) * -1;
    % Now double the angle, so that angles near pi and 0 get averaged
    % together correctly:
    dt = dt .* exp(1i * angle(dt));
    dt = dt(~isnan(dt));
    % Divide the angle by 2 to remove the doubling done on the previous
    % line.
    pang = angle(mean(dt,'all'))/2.0;
else
    pang = angle(mean(dt,'all'));
end

p_heading = round((90-rad2deg(pang)),4);

end

