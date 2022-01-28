function ds = set_declination(data_set, declin)
%   Set the magnetic declination
%
%    Parameters
%    ----------
%    declination : float
%       The value of the magnetic declination in degrees (positive
%       values specify that Magnetic North is clockwise from True North)
%
%    Returns
%    ----------
%    ds : dataset structure
%        Dataset adjusted for the magnetic declination
%
%    Notes
%    -----
%    This method modifies the data object in the following ways:
%
%    - If the dataset is in the *earth* reference frame at the time of
%      setting declination, it will be rotated into the "*True-East*,
%      *True-North*, Up" (hereafter, ETU) coordinate system
%
%    - dat['orientmat'] is modified to be an ETU to
%      instrument (XYZ) rotation matrix (rather than the magnetic-ENU to
%      XYZ rotation matrix). Therefore, all rotations to/from the 'earth'
%      frame will now be to/from this ETU coordinate system.
%
%    - The value of the specified declination will be stored in
%      dat.attrs['declination']
%
%    - dat['heading'] is adjusted for declination
%      (i.e., it is relative to True North).
%
%    - If dat.attrs['principal_heading'] is set, it is
%      adjusted to account for the orientation of the new 'True'
%      earth coordinate system (i.e., calling set_declination on a
%      data object in the principal coordinate system, then calling
%      dat.rotate2('earth') will yield a data object in the new
%      'True' earth coordinate system)
%  
    if isfield(data_set.attrs, 'declination')
        angle = declin - ds.attrs.declination;
    else
        angle = declin;
    end
    cd = cos(-deg2rad(angle));
    sd = sin(-deg2rad(angle));

    % The ordering is funny here because orientmat is the
    % transpose of the inst->earth rotation matrix:
    rdec = [cd, -sd, 0; sd, cd, 0; 0, 0, 1];

    if strcmpi(data_set.coord_sys,'earth')
        rotate2earth = true;
        data_set = rotate2(data_set, 'inst');
    else
        rotate2earth = false;
    end

    data_set.orientmat.data = tensorproduct(data_set.orientmat.data, ...
        rdec, 'edbc,ab->edac'); %'edcb,ab->edac' works for AWAC_test01.wpr
                                %'edbc,ab->edac' works for vector_data_imu01.VEC 

    if isfield(data_set,'heading')
        data_set.heading.data = data_set.heading.data + angle;
    end
    if rotate2earth
        data_set = rotate2(data_set, 'earth');
    end
    if isfield(data_set,'principal_heading')
        data_set.principal_heading = data_set.principal_heading + angle;
    end

    data_set.attrs.declination = declin;
    data_set.attrs.declination_in_orientmat = true;

    ds = data_set;
end

