function power_performance_wave(S, h, P, statistic, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     High-level function to compute the capture length matrix and Mean
%     Annual Energy Production (MAEP) for given wave spectra.
% 
% Parameters
% ------------
%   S: Spectral Density (m^2/Hz)
%       Pandas data frame
%           To make a pandas data frame from user supplied frequency and spectra
%           use py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(frequency,spectra)
%
%       OR
%
%       structure of form:
%           S.spectrum: Spectral Density (m^2/Hz)
%
%           S.type: String of the spectra type, i.e. Bretschneider, 
%           time series, date stamp etc.
%
%           S.frequency: frequency (Hz)
%
%   h: integer
%        Water depth [m]
%
%   P: array or vector
%        Power [W]
%
%   statistic: string or string array
%        Statistic for each bin, options include: "mean", "std", "median",
%        "count", "sum", "min", "max", and "frequency".  Note that "std"
%        uses a degree of freedom of 1 in accordance with IEC/TS 62600-100.
%        To output capture length matrices for multiple binning parameters,
%        define as a string array: statistic = ["", "", ""];
%        
%   rho: float (optional)
%        water density (kg/m^3)
%        to call: power_performance_wave(S,h,P,statistic,"rho",rho)
%
%   g: float (optional)
%        gravitational acceleration (m/s^2)
%        to call: power_performance_wave(S,h,P,statistic,"g",g)
%
%   frequency_bins: vector (optional) 
%      Bin widths for frequency of S. Required for unevenly sized bins
%
%    savepath: string (optional)
%         path and filename to save figure.
%         to call: power_performance_wave(S,h,P,statistic,"savepath",savepath)
%        
% Returns
% ---------
%   cl_matrix: figure
%       Capture length matrix
%
%   maep_matrix: float
%       Mean annual energy production
%
% Example
%     --------
%     >> S = read_NDBC_file('./data/wave/data.txt');
%     >> h = 60
%     >> P = randi([40,200], 743,1)
%     >> statistic = "mean"
%     >> power_performance_wave(S, h, P, statistic)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    S
    h
    P
    statistic
    options.rho = 1025;
    options.g = 9.80665;
    options.frequency_bins = "";
    options.savepath = "";
end

if ~isstruct(S)
    ME = MException('MATLAB:power_performance_wave','S must be a structure. See comments for formatting');
    throw(ME);
end

if any([~isnumeric(h), ~isnumeric(P)])
    ME = MException('MATLAB:power_performance_wave','h and P must be numbers');
    throw(ME);
end

% Compute the enegy periods from the NDBC spectra data 
Te = energy_period(S);

% Compute the significant wave height from the NDBC spectra data 
Hm0 = significant_wave_height(S);                                                

% Compute the energy flux  from the NDBC spectra data and water depth
J = energy_flux(S,h);

% calculating capture length with power and wave flux in vectors
L = capture_length(P,J);

% Need to set our Hm0 and Te bins for the capture length matrix
Hm0_bins = -0.5:0.5:max(fix(Hm0))+0.5; % Input is min, max, and n indecies for vector
Hm0_bins = Hm0_bins+0.25 ; 
Te_bins = 0:1:max(fix(Te));
Te_bins = Te_bins+0.5;

% Calculate the necessary capture length matrices for each statistic based
% on IEC/TS 62600-100
clmat.mean = capture_length_matrix(Hm0,Te,L,"mean",Hm0_bins,Te_bins);
clmat.std = capture_length_matrix(Hm0,Te,L,"std",Hm0_bins,Te_bins);
clmat.median = capture_length_matrix(Hm0,Te,L,"median",Hm0_bins,Te_bins);
clmat.count = capture_length_matrix(Hm0,Te,L,"count",Hm0_bins,Te_bins);
clmat.sum = capture_length_matrix(Hm0,Te,L,"sum",Hm0_bins,Te_bins);
clmat.min = capture_length_matrix(Hm0,Te,L,"min",Hm0_bins,Te_bins);
clmat.max = capture_length_matrix(Hm0,Te,L,"max",Hm0_bins,Te_bins);
clmat.freq = capture_length_matrix(Hm0,Te,L,"frequency",Hm0_bins,Te_bins);

% Capture Length Matrix using statistic
if length(statistic) > 1
    for i = 1:length(statistic)
        % Create wave energy flux matrix using statistic
        jmat(i) = wave_energy_flux_matrix(Hm0,Te,J,statistic(i),Hm0_bins,Te_bins);
        % Calcaulte MAEP from matrix
        maep_matrix(i) = mean_annual_energy_production_matrix(clmat.mean,jmat(i),clmat.freq);
        % Plot Capture Length Matrix
        if statistic(i) == "mean"
            clmat_out(i) = clmat.mean;
        elseif statistic(i) == "std"
            clmat_out(i) = clmat.std;
        elseif statistic(i) == "median"
            clmat_out(i) = clmat.median;
        elseif statistic(i) == "count"
            clmat_out(i) = clmat.count;
        elseif statistic(i) == "sum"
            clmat_out(i) = clmat.sum;
        elseif statistic(i) == "min"
            clmat_out(i) = clmat.min;
        elseif statistic(i) == "max"
            clmat_out(i) = clmat.max;
        elseif statistic(i) == "frequency"
            clmat_out(i) = clmat.freq;    
        else
            ME = MException('MATLAB:power_performance_wave','statistic must be a string or string array defined by one or multiple of the following: "mean", "std", "median","count", "sum", "min", "max", "frequency"');
            throw(ME);
        end
        figure(i)
        cl_matrix(i) = plot_matrix(clmat_out(i),"Capture Length");
        disp(cl_matrix(i))
    end
elseif length(statistic) == 1 
    if statistic == "mean"
        clmat = clmat.mean;
    elseif statistic == "std"
        clmat = clmat.std;
    elseif statistic == "median"
        clmat = clmat.median;
    elseif statistic == "count"
        clmat = clmat.count;
    elseif statistic == "sum"
        clmat = clmat.sum;
    elseif statistic == "min"
        clmat = clmat.min;
    elseif statistic == "max"
        clmat = clmat.max;
    elseif statistic == "frequency"
        clmat = clmat.frequency;
    else
        ME = MException('MATLAB:power_performance_wave','statistic must be a string or string array defined by one or multiple of the following: "mean", "std", "median","count", "sum", "min", "max", "frequency"');
        throw(ME);
    end
    cl_matrix = plot_matrix(clmat,"Capture Length");
    disp(cl_matrix)
end

len = strlength(options.savepath);
if len > 1
    saveas(figure, options.savepath);
end 

end

