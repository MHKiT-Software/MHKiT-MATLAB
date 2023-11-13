function figure=plot_flow_duration_curve(Q,F,options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Plots discharge vs exceedance probability as a Flow Duration Curve (FDC)
%
% Parameters
% ------------
%     Q: array
%         Discharge [m/s]
%
%     F: array
%          Exceedance probability [unitless]
%
%     title: string (optional)
%          title for the plot
%          to call: plot_flow_duration_curve(Q,F,"title",title)
%
%     savepath: string (optional)
%          path and filename to save figure.
%          to call: plot_flow_duration_curve(Q,F,"savepath",savepath)
%
% Returns
% ---------
%   figure: plot of discharge vs. exceedance probability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    Q
    F
    options.title = "";
    options.savepath = "";
end

temp.Q=Q;
temp.F=F;
T=struct2table(temp);
sortT=sortrows(T,'F','descend');


figure=plot(sortT.Q,sortT.F);
grid on
xlabel('Discharge [m^{3}/s]','FontSize',20)
ylabel('Exceedance Probability','FontSize',20)
set(gca,'XScale','log')

title(options.title)

len = strlength(options.savepath);
if len > 1
    saveas(figure, options.savepath);
end

hold off

