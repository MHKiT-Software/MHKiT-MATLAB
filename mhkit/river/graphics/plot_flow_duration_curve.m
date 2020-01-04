function figure=plot_flow_duration_curve(Q,F,varargin)
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
% 
% Returns
% ---------
%   figure: plot of discharge vs. exceedance probability 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

temp.Q=Q;
temp.F=F;
T=struct2table(temp);
sortT=sortrows(T,'F','descend');


figure=plot(sortT.Q,sortT.F);
grid on
xlabel('Discharge [m^{3}/s]','FontSize',20)
ylabel('Exceedance Probability','FontSize',20)
set(gca,'XScale','log')

if nargin == 3
    title(varargin{1})
end