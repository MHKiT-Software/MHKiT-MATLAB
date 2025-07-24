function DEL=damage_equivalent_load(data_signal,m,options)

%%%%%%%%%%%%%%%%%%%%
%     Calculates the damage equivalent load of a single data signal (or channel)
%     based on IEC TS 62600-3:2020 ED1. 4-point rainflow counting algorithm from
%     fatpack module is based on the following resources:
%
%     - `C. Amzallag et. al. Standardization of the rainflow counting method for
%       fatigue analysis. International Journal of Fatigue, 16 (1994) 287-293`
%       `ISO 12110-2, Metallic materials - Fatigue testing - Variable amplitude
%       fatigue testing.`
%     - `G. Marsh et. al. Review and application of Rainflow residue processing
%       techniques for accurate fatigue damage estimation. International Journal
%       of Fatigue, 82 (2016) 757-765`
%
% Parameters
% ------------
%     data_signal: vector
%         Data signal being analyzed
%
%     m : double or int
%         Fatigue slope factor of material
%
%     bin_num : int (optional)
%         Number of bins for rainflow counting method (minimum=100)
%         to call: get_DELs(data,chan_info,"bin_num",binNum)
%
%     data_length : double or int (optional)
%         Length of data in sec. Default for 1Hz is 600 seconds for 10min data
%         to call: get_DELs(data,chan_info,"data_length",t)
%
% Returns
% ---------
%     DEL: Structure
%         Damage equivalent load of signal
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    data_signal
    m
    options.bin_num {mustBeNumeric} = 100;
    options.data_length {mustBeNumeric} = 600;
end

data_signal_py = py.numpy.array(data_signal);
m_py = double(m);

DEL_py = py.mhkit.loads.general.damage_equivalent_load(data_signal_py, m_py, pyargs('bin_num',int32(options.bin_num),'data_length',int32(options.data_length)));

DEL = typecast_from_mhkit_python(DEL_py).data;
