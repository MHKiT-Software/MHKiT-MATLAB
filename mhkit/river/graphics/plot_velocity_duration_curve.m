function figure=plot_velocity_duration_curve(V,F,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Plots velocity vs exceedance probability as a Flow Duration Curve (FDC) 
%     
% Parameters
% ------------
%     V: array
%         Velocity [m/s] 
%         
%     F: array 
%          Exceedance probability [unitless]
%
%     title: string (optional)
%          title for the plot
%          to call: plot_velocity_duration_curve(P,F,"title",title)
%
%     savepath: string (optional)
%          path and filename to save figure.
%          to call: plot_velocity_duration_curve(P,F,"savepath",savepath)
% 
% Returns
% ---------
%   figure: plot of velocity vs. exceedance probability 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    V
    F
    options.title = "";
    options.savepath = "";
end

temp.V=V;

temp.F=F;

T=struct2table(temp);
sortT=sortrows(T,'F','descend');


figure=plot(sortT.V,sortT.F);
grid on
xlabel('Velocity [m/s]','FontSize',20)
ylabel('Exceedance Probability','FontSize',20)


title(options.title)

len = strlength(options.savepath);
if len > 1
    saveas(figure, options.savepath);
end 

hold off