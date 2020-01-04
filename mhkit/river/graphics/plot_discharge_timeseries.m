function figure=plot_discharge_timeseries(Q,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Plots discharge vs time 
%     
% Parameters
% ------------
%     Q: structure
%
%      Q.Discharge=Discharge [m/s]
%
%      Q.time= time
%         
%     title: string (optional)
%       title for the plot
% 
% Returns
% ---------
%     figure: Plot of discharge vs. time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure=plot(Q.time,Q.Discharge);
grid on
xlabel('Time','FontSize',20)
ylabel('Discharge [m^{3}/s]','FontSize',20)

if nargin ==2
    title(varargin{1})
end