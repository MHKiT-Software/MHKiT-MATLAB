function figure=plot_spectrum(wave_spectra)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Plots wave amplitude spectrum
%     
% Parameters
% ----------
%     wave_spectra: Structure of the following form:
%
%         wave_spectra.spectrum: Spectral Density (m^2-s;
%
%         wave_spectra.type: String of the spectra type, i.e. Bretschneider, time series, date stamp etc.
%
%         wave_spectra.frequency: frequency (Hz);
%         
% Returns
% ---------
%     figure: figure
%         Plot of wave amplitude spectra versus omega
%     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure=plot(wave_spectra.frequency*2*3.14,wave_spectra.spectrum/(2*3.14));
xlabel('Omega (^{rad}/_{s})')
ylabel('Spectral Density (m^{2}s/_{rad})') 
si=size(wave_spectra.spectrum);
if si(1)==1
    Hm0=significant_wave_height(wave_spectra);
    format_Spec='Spectrum: %s, Tp= %f, Hm0= %f';
    title(sprintf(format_Spec,wave_spectra.type,wave_spectra.Tp,Hm0));
end
    