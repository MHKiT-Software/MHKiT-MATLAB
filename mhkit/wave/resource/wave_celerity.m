function Cg=wave_celerity(k,h,options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates wave celerity (group velocity)
%
% Parameters
% ------------
%    k: wave number (1/m)
%        structure of form:
%           k.values= wave number
%
%           k.frequency= frequency (Hz)
%
%    h: double
%         Water depth (m)
%
%    g: double (optional)
%         gravitational acceleration (m/s^2)
%         to call: energy_flux(k,h,"g",g)
%
%    depth_check: bool (optional)
%         If True check depth regime. Default False.
%         to call: energy_flux(k,h,"depth_check",py.True)
%
%    ratio: double or int (optional)
%         Only applied if depth_check=True. If h/l > ratio,
%         water depth will be set to deep. Default ratio = 2
%         to call: energy_flux(k,h,"ratio",ratio)
%
% Returns
% -------
%     Cg: structure
%
%
%       Cg.values: water celerity
%
%       Cg.frequency [Hz]
%
%       Cg.h: height [m]
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    k
    h
    options.g = 9.80665;
    options.depth_check = py.False;
    options.ratio = 2;

end

if (isa(k,'py.pandas.core.frame.DataFrame')~=1)
    if (isstruct(k)==1)

            k=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(k.frequency,py.numpy.array(k.values),1);

    else
        ME = MException('MATLAB:wave_celerity','k needs to be a structure or Pandas dataframe, use py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas to create one');
        throw(ME);
    end
end


Cgdf=py.mhkit.wave.resource.wave_celerity(k,h,pyargs('g',options.g,...
    'depth_check',options.depth_check,'ratio',options.ratio));

Cg.values=double(py.array.array('d',py.numpy.nditer(Cgdf.values)));
Cg.frequency=double(py.array.array('d',py.numpy.nditer(Cgdf.index)));
Cg.h=h;

