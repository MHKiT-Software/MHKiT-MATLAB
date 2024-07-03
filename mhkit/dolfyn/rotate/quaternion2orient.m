function omat = quaternion2orient(quaternions)
%     Calculate orientation from Nortek AHRS quaternions, where
%       q = [W, X, Y, Z] instead of the standard
%       q = [X, Y, Z, W] = [q1, q2, q3, q4]
%
%     Parameters
%     ----------
%     quaternions : Structure
%         Quaternion structure from the raw dataset
%
%     Returns
%     -------
%     orientmat : |ndarray|
%         The inst2earth rotation maxtrix as calculated from the
%         quaternions
%
%     See Also
%     --------
%     scipy.spatial.transform.Rotation

    omat = struct();
    omat.data = zeros([length(quaternions.coords.time),1,3,3], 'double');
    omat.dims = { 'time', 'inst', 'earth' };
    omat.coords.time = quaternions.coords.time;
    omat.coords.inst = {'X' 'Y' 'Z'};
    omat.coords.earth = {'E' 'N' 'U'};

    for i = 1:length(quaternions.coords.time)
        r = quat2rotmat([quaternions.data(i,:,2),...
                         quaternions.data(i,:,3),...
                         quaternions.data(i,:,4),...
                         quaternions.data(i,:,1)]);
        omat.data(i,:,:,:) = r;
    end

    omat.data = permute(omat.data,[1, 2, 4, 3]);

    function R = quat2rotmat(quat)
        R = zeros([3,3]);
        quat = normalize(quat,'norm');

        R(1,1) =  1-2*(quat(2)^2 + quat(3)^2);
        R(2,1) = -2*(quat(3)*quat(4) - quat(2)*quat(1));
        R(3,1) =  2*(quat(2)*quat(4) + quat(3)*quat(1));

        R(1,2) =  2*(quat(3)*quat(4) + quat(2)*quat(1));
        R(2,2) = -1*(1-2*(quat(2)^2 + quat(4)^2));
        R(3,2) =  2*(quat(2)*quat(3) - quat(4)*quat(1));

        R(1,3) = -2*(quat(2)*quat(4) - quat(3)*quat(1));
        R(2,3) =  2*(quat(2)*quat(3) + quat(4)*quat(1));
        R(3,3) = -1*(1-2*(quat(3)^2 + quat(4)^2));
    end

end

