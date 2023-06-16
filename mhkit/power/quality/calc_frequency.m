function [alpha0,freq] = calc_frequency(u_m,methodopts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculates the fundamental frequency of the measured voltage u_m(t)
%   using Short-Time Fourier Transform (STFT).
%   
% Parameters
% -----------
%   u_m: struct()
%       .time: time for each time step; 
%       .data: measured voltage (V) with instantaneous values u_m(t).
%   methodopts: cell array 
%       cell array of name-value arguments, e.g.,
%       opts = {'Window',rectwin(M),'OverlapLength',L,'FFTLength',128}
%       
% 
% Returns
% -------
%   freq: struct()
%       .time: time instants corresponding to the centers of the data 
%              segments used to capture the power spectrum estimates. 
%       .data: freq(t), the fundamental frequency calculated from the 
%              STFT for u_m.
%   alpha0: double 
%       electrical angle at t=0.
% Note
% -------
% 1. Short-Time Fourier Transform (STFT) is used to calculate the frequency.
% calculate frequency using u_m and corresponding method:
%    [alpha0,freq] = calc_frequency(u_m,methodopts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    freq = struct();
    d = seconds(u_m.time(2)-u_m.time(1));
    % what is window size one cycle?
    if isempty(methodopts)
        fprintf('calc_frequency: using default settings of stft.\n');
        [s,f,t]=stft(u_m.data,d);
    else
        [s,f,t]=stft(u_m.data,d,methodopts{:});
    end   
    m = abs(s);
    [~,fidx] = max(m,[],1);
    f_fund = f(fidx);
    freq.time = t;
    freq.data = f_fund;
    fftout = fft(u_m.data);
    alpha0 = angle(fftout(1));
end