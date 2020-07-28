function figure=plot_power_duration_curve(P,F,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Plots power vs exceedance probability as a Flow Duration Curve (FDC) 
%     
% Parameters
% ------------
%     P: array
%         Power [W] 
%         
%     F: array 
%          Exceedance probability [unitless]
%
%     title: string (Optional)
%          title for the plot
%          to call: plot_power_duration_curve(P,F,"title",title)
%
%     savepath: string (optional)
%          path and filename to save figure.
%          to call: plot_power_duration_curve(P,F,"savepath",savepath)
% 
% Returns
% ---------
%   figure: plot of power duration curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    P
    F
    options.title = "";
    options.savepath = "";
end

temp.P=P;

temp.F=F;

T=struct2table(temp);
sortT=sortrows(T,'F','descend');


figure=plot(sortT.P,sortT.F);
grid on
xlabel('Power [W]','FontSize',20)
ylabel('Exceedance Probability','FontSize',20)


title(options.title)

len = strlength(options.savepath);
if len > 1
    saveas(figure, options.savepath);
end 

hold off