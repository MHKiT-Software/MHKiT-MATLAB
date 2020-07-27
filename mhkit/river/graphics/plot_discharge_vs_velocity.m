function figure=plot_discharge_vs_velocity(Q,V,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Plots discharge vs velocity 
%     
% Parameters
% ------------
%     Q: array
%         Discharge [m/s]
%
%     V: array
%         Velocity [m/s] 
%
%     title: string (optional)
%         title for the plot
%         to call: plot_discharge_vs_velocity(Q,V,"title",title)
%
%     savepath: string (optional)
%         path and filename to save figure.
%         to call: plot_discharge_vs_velocity(Q,V,"savepath",savepath)
%
%     polynomial_coeff: array (optional)
%         polynomial coefficients which can be computed from 
%         polynomial_fit.m. Expects poly.coef
%         to call: plot_discharge_vs_velocity(Q,V,"polynomial_coeff",polynomial_coeff)
% 
% Returns
% ---------
%   figure: plot of discharge vs. velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    Q
    V
    options.title = "";
    options.savepath = "";
    options.polynomial_coeff = [];
end

figure=plot(Q,V,'o');
hold on
grid on
xlabel('Discharge [m^{3}/s]','FontSize',20)
ylabel('Velocity [m/s]','FontSize',20)

title(options.title)
mi=min(Q);
ma=max(Q);
len = strlength(options.savepath);
lenp = length(options.polynomial_coeff);
if lenp > 1
   x=linspace(mi(1),ma(1),100);
   p2=plot(x,polyval(options.polynomial_coeff,x),'--');
   legend(p2,{'Polynomial Fit'},'Location','northwest') 
end

if len > 1
    saveas(figure, options.savepath);
end 



hold off      
        
