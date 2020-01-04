function figure=plot_power_duration_curve(P,F,varargin)
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
%       title for the plot
% 
% Returns
% ---------
%   figure: plot of power vs. exceedance probability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

temp.P=P;

temp.F=F;

T=struct2table(temp);
sortT=sortrows(T,'F','descend');


figure=plot(sortT.P,sortT.F);
grid on
xlabel('Power [W]','FontSize',20)
ylabel('Exceedance Probability','FontSize',20)


if nargin == 3
    title(varargin{1})
end