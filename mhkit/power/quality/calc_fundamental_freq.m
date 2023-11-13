function [alpha0,freq] = calc_fundamental_freq(u_m,method,methodopts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculates the fundamental frequency of the measured voltage u_m(t)
%   using Zero-Crossing-Detection (ZCD) or Short-Time Fourier
%   Transform (STFT).
%
% Parameters
% -----------
%   u_m: struct()
%       .time: time for each time step;
%       .data: array of size (ntime), measured voltage (V) with instantaneous values u_m(t).
%
%   method: string (case-insensitive)
%       can be either 'zcd' or 'stft'.
%       'zcd': Zero-Crossing-Detection method
%       'stft': Short-Time Fourier Transform method
%
%   methodopts: cell array
%       Name-value arguments for STFT, e.g.,
%       opts = {'Window',rectwin(M),'OverlapLength',L,
%               'FFTLength',Nfrq,'FrequencyRange','onesided'}
%       if using 'ZCD' method, can set to {}.
% Returns
% -------
%   freq: struct()
%       .time: time at each measurement time step (s)
%       .data: double array of size (ntime), freq(t), the fundamental
%       frequency (that may vary over time) calculated for u_m(t).
%   alpha0: double
%       The electrical angle (radians) at t=0.
% Note
% -------
% 1. Window settings are crucial to get accurate results if
% Short-Time Fourier Transform (STFT) is used to calculate the frequency.
% 2. Steps:
%   Step1. adjust time values
%   Step2. Estimate fundamental freq using either STFT or ZCD
%   Step3. Interpolate estimated freq back to u_m.time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % check input:
    if ~isfield(u_m,'time') || ~isfield(u_m, 'data')
        ME = MException('MATLAB:calc_fundamental_freq',...
            'invalid handles in structure, must contain x.data & x.time');
        throw(ME);
    end
    % 1. make sure time starts at 0.0
    u_m.time = u_m.time-u_m.time(1);
    freq = struct();
    d = seconds(u_m.time(2)-u_m.time(1));
    % 2. Estimate fundamental frequency using STFT or ZCD
    % opt 1, use STFT
    if strcmpi(method,'stft')
        % what is window size one cycle?
        if isempty(methodopts)
            fprintf('calc_frequency: using default settings of stft.\n');
            [s,f,t]=stft(u_m.data,d);
        else
            [s,f,t]=stft(u_m.data,d,methodopts{:});
        end
        m = abs(s);
        [~,fidx] = max(m,[],1);
        f_mid = f(fidx);
        t_lead = seconds([2*t(1)-t(2);0.5*(t(1:end-1)+t(2:end))]);%t;
        amp = min(abs(max(u_m.data)),abs(min(u_m.data)));
        if u_m.data(1)==0 && u_m.data(2)>0
            alpha0 = 0;
        elseif u_m.data(1)==0 && u_m.data(2)<0
            alpha0 = pi;
        else
            alpha0 = asin(u_m.data(1)/amp);
        end
    else
        % opt 2, use ZCD
        x = u_m.time; y = u_m.data;
        % indexes of zero crossings: x(i)*x(i+1)<0
        i=find(y(1:end-1).*y(2:end)<0);
        m = (y(i+1)-y(i))./(x(i+1)-x(i));
        x0 = -y(i)./m+x(i);
        t_lead = x0(1:end-1);%0.5*(x0(1:end-1)+x0(2:end));
        f_mid = 1./(diff(x0)*2);
        if u_m.data(1)==0 && u_m.data(2)>0
            alpha0 = 0;
        elseif u_m.data(1)==0 && u_m.data(2)<0
            alpha0 = pi;
        else
            alpha0 = pi*(1+sign(u_m.data(1)))/2-2*pi*f_mid(1)*x0(1);
        end
    end
    % 3. Interpolate estimated frequency back to u_m.time resolution
    ft = interp1(t_lead,f_mid,u_m.time,'previous','extrap');
    ft(isnan(ft))=f_mid(1);
    freq.data = ft;
    freq.time = u_m.time;

end

