function figure=plot_matrix(M,Mtype, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Plots the matrix with Hm0 and Te on the y and x axis
%
% Parameters
% ----------
%    M: structure
%
%         M.values: matrix
%
%         M.Hm0_bins
%
%         M.Te_bins
%
%         M.stat
%
%    Mtype: string
%         type of matrix (i.e. power, capture length, etc.) to be used
%         in plot title 
%
%     savepath: string (optional)
%         path and filename to save figure.
%         to call: plot_matrix(M, Mtype,"savepath",savepath)
%
% Returns
% ---------
%   figure: plot of the matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    M 
    Mtype
    options.savepath = "";
end

figure=pcolor(M.Te_bins,M.Hm0_bins,M.values);
colormap(flipud(hot))
ylabel('Hm0 [m]','FontSize',20)
xlabel('Te [s]','FontSize',20)
x=strcat(Mtype,': ',M.stat);
title(x)
colorbar
pos=get(gca,'position');
[rows,cols]=size(M.values);
width=pos(3)/(cols-1);
height =pos(4)/(rows-1);
%create textbox annotations
% for i=1:cols-1
%       for j=rows-1:-1:1 
%           if ~isnan(M.values(j,i)) 
%         annotation('textbox',[pos(1)+width*(i-1),pos(2)+height*(j-1),width,height], ...
%        'string',num2str(M.values(j,i)),'LineStyle','none','HorizontalAlignment','center',...
%        'VerticalAlignment','middle');
%           end
%       end
% end

len = strlength(options.savepath);
if len > 1
    saveas(figure, options.savepath);
end 
