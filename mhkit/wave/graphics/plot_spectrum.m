function figure=plot_spectrum(wave_spectra,options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Plots wave amplitude spectrum
%
% Parameters
% ----------
%	wave_spectra: Structure of the following form:
%
%       wave_spectra.spectrum: Spectral Density (m^2-s;
%
%       wave_spectra.type: String of the spectra type, i.e. Bretschneider, time series, date stamp etc.
%
%       wave_spectra.frequency: frequency (Hz);
%
%     savepath: string (optional)
%         path and filename to save figure.
%         to call: plot_spectrum(wave_spectra,"savepath",savepath)
%
% Returns
% ---------
%	figure: figure
%       Plot of wave amplitude spectra versus omega
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    wave_spectra
    options.savepath = "";
end

figure=plot(wave_spectra.frequency*2*3.14,wave_spectra.spectrum/(2*3.14));
xlabel('Omega (^{rad}/_{s})')
ylabel('Spectral Density (m^{2}s/_{rad})')
si=size(wave_spectra.spectrum);
if si(1)==1
    Hm0=significant_wave_height(wave_spectra);
    Tp=peak_period(wave_spectra);
    format_Spec='Spectrum: %s, Tp= %f, Hm0= %f';
    title(sprintf(format_Spec,wave_spectra.type,Tp,Hm0));
end

len = strlength(options.savepath);
if len > 1
    saveas(figure, options.savepath);
end
