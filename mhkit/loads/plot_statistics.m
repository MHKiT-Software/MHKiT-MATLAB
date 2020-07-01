function figure=plot_statistics(x,y_mean,y_min,y_max,options)

%%%%%%%%%%%%%%%%%%%%
%     Plot showing standard raw statistics of variable 
%     
% Parameters
% ------------
%     x: vector  
%         vector of x-axis values
%
%     y_mean: vector
%         Vector of mean values
%
%     y_max: vector
%         Vector of max values
%
%     y_min: vector
%         Vector of min values
%
%     y_stdev: vector (optional)
%         Vector of stand drviations 
%         to call: statplotter(x,y_mean,y_max,y_min,"y_stdev",vstdev)
%
%     xlabel: string (optional)
%         x-axis lable for the plot
%         to call: statplotter(x,y_mean,y_max,y_min,"xlable",xlable)
%
%     ylabel: string (optional)
%         y-axis lable for the plot
%         to call: statplotter(x,y_mean,y_max,y_min,"ylable",ylable)
%
%     title: string (optional)
%         title for the plot
%         to call: statplotter(x,y_mean,y_max,y_min,"title",title)
%
%     savepath: string (optional)
%         path and filename to save figure.
%         to call: statplotter(x,y_mean,y_max,y_min,"savepath",savepath)
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
    y_mean
    y_min
    y_max
    options.y_stdev = [];
    options.xlabel = "";
    options.ylabel = "";
    options.title = "";
    options.savepath = "";
end

figure = scatter(x,y_mean);
hold on
grid on
scatter(x,y_max)
scatter(x,y_min)
if size(options.y_stdev) == size(x)
    scatter(x,options.y_stdev)
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


