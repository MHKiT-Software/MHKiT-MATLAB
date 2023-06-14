function [alpha0,freq] = calc_frequency(u_m,methodopts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculates the frequency of the measured voltage u_m(t).
%   
% Parameters
% -----------
%   u_m: struct()
%       .time: time for each time step; 
%       .data: measured voltage with instantaneous value u_m(t) (V).
%   methodopts: cell array 
%       cell array of STFT options, e.g.,
%       opts = {'Window',rectwin(10),'OverlapLength',8,'FFTLength',128}
%       
% 
% Returns
% -------
%   freq: struct()
%       .time: time for each time step; 
%       .data: freq(t), the frequency for measured voltage.
% Note
% -------
%   Short-Time Fourier Transform (STFT) is used to calculate the frequency.
%   
%   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate frequency using STFT or ZCD
    
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
    f_fund = zeros(length(u_m.data),1);
    f_fund(f(fidx);
    freq.time = t;
    freq.data = f_fund;
    % u_m.time(1) may !=0:
    % alpha0 = angle(s(fidx(1),1))-2*pi*f_fund(1)*seconds(u_m.time(1)-t(1));
    fftout = fft(u_m.data);
    alpha0 = angle(fftout(1));
end