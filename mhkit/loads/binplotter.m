function figure=binplotter(bcenters,bmean,bmax,bmin,bstdmean,bstdmax,bstdmin,options)

%%%%%%%%%%%%%%%%%%%%
%     plot showing standard binned statistics of single variable 
%     
% Parameters
% ------------
%     bcenters: vector  
%         x-axis bin center values
%
%     bmean: vector
%         binned mean statistical values of variable
%
%     bmax: vector
%         binned max statistical values of variable
%
%     bmin: vector
%         binned min statistical values of variable
%
%     bstdmean: vector 
%         standard deviations of mean binned statistics
%
%     bstdmax: vector 
%         standard deviations of max binned statistics
%
%     bstdmin: vector 
%         standard deviations of min binned statistics
%
%     xlabel: string (optional)
%         x-axis lable for the plot
%         to call: binplotter(bcenters,bmean,bmax,bmin,bstdmean,bstdmax,bstdmin,"xlable",xlable)
%
%     ylabel: string (optional)
%         y-axis lable for the plot
%         to call: binplotter(bcenters,bmean,bmax,bmin,bstdmean,bstdmax,bstdmin,"ylable",ylable)
%
%     title: string (optional)
%         title for the plot
%         to call: binplotter(bcenters,bmean,bmax,bmin,bstdmean,bstdmax,bstdmin,"title",title)
%
%     savepath: string (optional)
%         path and filename to save figure.
%         to call: binplotter(bcenters,bmean,bmax,bmin,bstdmean,bstdmax,bstdmin,"savepath",savepath)
%
% Returns
% ---------
%     figure: figure 
%         scatter plot of statistics 
%
%         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    bcenters 
    bmean
    bmax
    bmin
    bstdmean
    bstdmax
    bstdmin
    options.xlabel = "";
    options.ylabel = "";
    options.title = "";
    options.savepath = "";
end
I = ~isnan(bmax);
figure = errorbar(bcenters,bmean,bstdmean);
hold on
grid on
errorbar(bcenters(I),bmax(I),bstdmax(I));
errorbar(bcenters(I),bmin(I),bstdmin(I));
legend({'mean','max','min'},'Location','northeast')


title(options.title)
xlabel(options.xlabel)
ylabel(options.ylabel)
hold off