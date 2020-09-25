function figure=plot_environmental_contours(x1, x2, x1_contour, x2_contour, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Plots an overlay of the x1 and x2 variables to the calculated
%     environmental contours.
%     
% Parameters
% ------------
%     x1: vector
%         x-axis data
%
%     x2: vector
%         y-axis data
%
%     x1_countour: Table
%         Calculated x1 contour values
%
%     x2_countour: Table
%         Calculated x2 contour values
%
%     x_label: string (optional)
%         x-axis label. Default None.
%         to call: plot_environmental_contours(x1,x2,x1_contour,x2_contour,"x_label",x_label)
%
%     y_label: string (optional)
%         y-axis label. Default None.
%         to call: plot_environmental_contours(x1,x2,x1_contour,x2_contour,"y_label",y_label)
%
%     data_label: string (optional)
%         Legend label for x1, x2 data (e.g. 'Buoy 46022'). 
%         Default None.
%         to call: plot_environmental_contours(x1,x2,x1_contour,x2_contour,"data_label",data_label)
%
%     contour_label: string or array of strings (optional)
%         Legend label for x1_contour, x2_contour countor data 
%         (e.g. '100-year contour'). Default None.
%         to call: plot_environmental_contours(x1,x2,x1_contour,x2_contour,"contour_label",contour_label)
%
%     title: string (optional)
%         title for the plot
%         to call: plot_environmental_contours(x1,x2,x1_contour,x2_contour,"title",title)
%
%     savepath: string (optional)
%         path and filename to save figure.
%         to call: plot_environmental_contours(x1,x2,x1_contour,x2_contour,"savepath",savepath)
%         
% Returns
% ---------
%     figure: figure
%         Envoronmental contour plot
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    x1
    x2
    x1_contour
    x2_contour
    options.x_label = "";
    options.y_label = "";
    options.data_label = "";
    options.contour_label ="";
    options.title = "";
    options.savepath = "";
end

figure=scatter(x1,x2);
hold on
for field = fieldnames(x1_contour)'
    if ~any(strcmp(field{1} , {'Properties','Row','Variables'} ))
    p2 = plot(x1_contour.(field{1}), x2_contour.(field{1}));
    end
end
grid on
xlabel(options.x_label)
ylabel(options.y_label) 

title(options.title)

len = strlength(options.savepath);
if len > 1
    saveas(figure, options.savepath);
end 

hold off