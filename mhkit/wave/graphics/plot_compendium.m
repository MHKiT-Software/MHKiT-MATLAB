function f = plot_compendium(Hs, Tp, Dp, time, options)
%PLOT_COMPENDIUM Creates subplots of wave height, peak period and direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Creates subplots showing: Significant Wave Height (Hs),
%   Peak Period (Tp), and Direction (Dp). Developed based on:
%   http://cdip.ucsd.edu/themes/media/docs/documents/html_pages/compendium.html
%   
%   Parameters
%   ----------
%       Hs: double array
%           significant wave height
%       Tp: double array
%           peak period
%       Dp: double array
%           direction
%       time: double array
%           timestamps as POSIX
%       buoy_title: string (optional)
%           figure suptitle (super title)
%   
%   Returns
%   -------
%       f : figure object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    Hs (1,:) double
    Tp (1,:) double
    Dp (1,:) double
    time (1,:) double
    options.buoy_title string = "";
end

timestamps = datetime(time, 'convertfrom','posixtime');
title_text = sprintf("%s to %s", ...
    datestr(timestamps(1), 'yyyy-mm-dd'), ...
    datestr(timestamps(end), 'yyyy-mm-dd'));

% Create figure window
f = figure( ...
    'Name', 'Compendium', ...
    'visible', 'off');                          % hide before resizing
f.Position(3:4) = 2 * f.Position(3:4);          % double default size
movegui(f, 'center');
set(f, 'visible', 'on');

% Hs subplot
ax1 = subplot(3, 1, 1);
yyaxis left
plot(timestamps, Hs)
title(title_text, "FontSize", 15)
grid on
ylabel('Hs [m]', "FontSize", 14)
set(ax1, 'xticklabel', [], 'Ycolor', 'k')
ylim_left = ylim;
yyaxis right                                    % add right y-axis
ft_in_meters = 3.28084;
set(ax1, 'ylim', ft_in_meters * ylim_left, 'Ycolor', 'k')
ylabel('Hs [ft]', "FontSize", 14)

% Tp subplot
ax2 = subplot(3, 1, 2);
plot(timestamps, Tp)
grid on
set(ax2, 'xticklabel', [])                      % hide x-axis labels
ylabel('Tp [s]', "FontSize", 14)

% Dp subplot
ax3 = subplot(3, 1, 3);
scatter(timestamps, Dp, 10, "filled")
ylim([0 360])
yticks(0:90:360)
grid on
xlabel('Date', "FontSize", 14)
ylabel('Dp [deg]', "FontSize", 14)

% Add super title over all sublots
if options.buoy_title ~= ""
    sgtitle('Test Super Title', "FontSize", 20)
end
