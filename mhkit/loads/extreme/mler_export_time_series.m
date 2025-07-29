function mler_ts = mler_export_time_series(RAO, mler, sim, k)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Generate the wave amplitude time series at X0 from the calculated
%     MLER coefficients
%
%     Parameters
%     ----------
%     RAO: array
%         Response amplitude operator.
%     mler: struct
%         MLER coefficients dataframe generated from an MLER function.
%     sim: struct
%         Simulation parameters formatted by output from
%         'mler_simulation'.
%     k: array
%         Wave number.
%
%     Returns
%     -------
%     mler_ts: struct
%         Contains resulting wave height [m] and linear response [*] indexed
%         by time [s].
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assert statements
if ~isa(RAO,'numeric')
    error('ERROR: RAO must be an array')
end
if ~isa(mler,'struct')
    error('ERROR: mler must be a struct')
end
if ~isa(sim,'struct')
    error('ERROR: sim must be a struct')
end
if ~isa(k,'numeric')
    error('ERROR: k must be an array')
end

freq = mler.frequency * 2*pi;  % convert Hz to rad/s
dw = (max(freq) - min(freq)) / (length(freq)-1);  % get delta

% calculate the series
wave_amp_time = zeros(sim.maxIT, 2);
xi = sim.X0;
for i=1:length(sim.T)
    % conditioned wave
    wave_amp_time(i,1) = sum( ...
        sqrt(mler.conditioned_spectrum.*2*dw) .* ...
        cos(freq.*(sim.T(i)-sim.T0) + mler.phase - k'.*(xi-sim.X0)));
    % Response calculation
    wave_amp_time(i,2) = sum( ...
        sqrt(mler.conditioned_spectrum.*2*dw) .* abs(RAO) .* ...
        cos(freq.*(sim.T(i)-sim.T0) - k'.*(xi-sim.X0)));

% generate outputs
mler_ts.wave_height = wave_amp_time(:,1);
mler_ts.linear_response = wave_amp_time(:,2);
mler_ts.time = sim.T;

end
