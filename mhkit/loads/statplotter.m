function figure=statplotter(x,vmean,vmin,vmax,options)

%%%%%%%%%%%%%%%%%%%%
%     plot showing standard raw statistics of variable 
%     
% Parameters
% ------------
%     x: vector  
%         vector of x-axis values
%
%     vmean: vector
%         vector of mean values
%
%     vmax: vector
%         vector of max values
%
%     vmin: vector
%         vector of min values
%
%     vstdev: vector (optional)
%         vector of stand drviations 
%         to call: statplotter(x,vmean,vmax,vmin,"vstdev",vstdev)
%
%     xlabel: string (optional)
%         x-axis lable for the plot
%         to call: statplotter(x,vmean,vmax,vmin,"xlable",xlable)
%
%     ylabel: string (optional)
%         y-axis lable for the plot
%         to call: statplotter(x,vmean,vmax,vmin,"ylable",ylable)
%
%     title: string (optional)
%         title for the plot
%         to call: statplotter(x,vmean,vmax,vmin,"title",title)
%
%     savepath: string (optional)
%         path and filename to save figure.
%         to call: statplotter(x,vmean,vmax,vmin,"savepath",savepath)
%
% Returns
% ---------
%     figure: figure 
%         scatter plot of statistics 
%
%         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    x 
    vmean
    vmin
    vmax
    options.vstdev = [];
    options.xlabel = "";
    options.ylabel = "";
    options.title = "";
    options.savepath = "";
end

figure = scatter(x,vmean);
hold on
grid on
scatter(x,vmax)
scatter(x,vmin)
if size(options.vstdev) == size(x)
    scatter(x,options.vstdev)
    legend({'mean','max','min','stdev'},'Location','northeast')
else
    legend({'mean','max','min'},'Location','northeast')
end

title(options.title)
xlabel(options.xlabel)
ylabel(options.ylabel)

len = strlength(options.savepath);
if len > 1
    saveas(figure, options.savepath);
end 

hold off


