function figure=plot_discharge_vs_velocity(Q,V,varargin)
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
%
%   polynomial_coeff: array (optional)
%       polynomial coefficients which can be computed from 
%       polynomial_fit.m. Expects poly.coef
% 
% Returns
% ---------
%   figure: plot of discharge vs. velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure=plot(Q,V,'o');
grid on
xlabel('Discharge [m^{3}/s]','FontSize',20)
ylabel('Velocity [m/s]','FontSize',20)

if nargin == 3
    if isstring(varargin{1})
        title(varargin{1})
    elseif isnumeric(varargin{1})
        mi=min(Q);
        
        ma=max(Q);
        si=size(varargin{1});
        x=linspace(mi(1),ma(1),100);
        hold on
        p2=plot(x,polyval(varargin{1},x),'--');
        legend(p2,{'Polynomial Fit'},'Location','northwest')
    else
        ME = MException('MATLAB:plot_discharge_vs_velocity','Variable argument is of the wrong type');
        throw(ME);
    end
    
elseif nargin ==4 
    mi=min(Q);
    ma=max(Q);
    if isstring(varargin{1}) & isnumeric(varargin{2})
        
        si=size(varargin{2});
        x=linspace(mi(1),ma(1),100);
        hold on
        p2=plot(x,polyval(varargin{2},x),'--');
        legend(p2,{'Polynomial Fit'},'Location','northwest')
        title(varargin{1})
    elseif isstring(varargin{2}) & isnumeric(varargin{1})
        si=size(varargin{1});
        x=linspace(mi(1),ma(1),100);
        hold on
        p2=plot(x,polyval(varargin{1},x),'--');
        legend(p2,{'Polynomial Fit'},'Location','northwest')
        title(varargin{2})
    else
        ME = MException('MATLAB:plot_discharge_vs_velocity','Variable argument is of the wrong type');
        throw(ME);
    end
elseif nargin > 4
    ME = MException('MATLAB:plot_discharge_vs_velocity','Too many arguments given');
        throw(ME);
end
        
hold off        
        
