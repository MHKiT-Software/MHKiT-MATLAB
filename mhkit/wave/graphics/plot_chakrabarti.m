function figure=plot_chakrabarti(H, lambda_w, D, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Plots, in the style of Chakrabarti (2005), relative importance of viscous,
%     inertia, and diffraction phemonena
% 
%     Chakrabarti, Subrata. Handbook of Offshore Engineering (2-volume set).
%     Elsevier, 2005.
% 
% Parameters
% ------------
%    H: integer, double or vector
%         Wave height [m]
%
%    lambda_w: integer, double or vector
%         Wave length [m]
%
%    D: integer, double or vector of 
%         Characteristic length [m]
%
%    savepath: string (optional)
%         path and filename to save figure.
%         to call: plot_chakrabarti(H,lambda_w,D,"savepath",savepath)
%         
% Returns
% ---------
%	 figure: figure
%         Plots wave force regime as Keulegan-Carpenter parameter versus
%         diffraction parameter
%
% Examples
%     --------
%     **Using Integers**
%     >> D = 5
%     >> H = 8
%     >> lambda_w = 200
%     >> plot_chakrabarti(H,lambda_w,D)
% 
%     **Using vector**
%     >> D = linspace(5,15,5)
%     >> H = 8*ones(size(D))
%     >> lambda_w = 200*ones(size(D))
%     >> plot_chakrabarti(H,lambda_w,D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    H
    lambda_w
    D
    options.savepath = "";
end

if any([~isnumeric(H), ~isnumeric(lambda_w), ~isnumeric(D)])
    ME = MException('MATLAB:plot_chakrabarti','H,lambda_w, and D must be numbers');
    throw(ME);
end

if any([length(H) ~= length(lambda_w), length(lambda_w) ~= length(D)])
    ME = MException('MATLAB:plot_chakrabarti','H,lambda_w, and D must all be the same length');
    throw(ME);
end

szD = size(D);
KC = zeros(szD);
Diffraction = zeros(szD);
for i=1:max(szD)
    KC(i) = H(i) ./ D(i);
    Diffraction(i) = pi.*D(i) ./ lambda_w(i);
    txt = ('H='+string(H(i))+', {\lambda_w}='+string(lambda_w(i))+', D='+string(D(i)));
    figure = loglog(Diffraction(i), KC(i),'o', 'DisplayName', txt);
    legend('Location', 'southoutside', 'NumColumns', 2);
    xlabel('Diffraction parameter, $\frac{\pi*D}{\lambda_w}$','Interpreter','latex')
    ylabel('KC parameter, $\frac{H}{D}$','Interpreter','latex')
    set(figure, 'markerfacecolor', get(figure, 'color'))
    hold on
end

if any(KC>=10) || any(KC<=.02) || any(Diffraction>=50) || any(lambda_w >= 1000)
    axis auto

else
    xlim([0.01 10])
    ylim([0.01 50])
    
end
xDrag = 0.0125;
xInertiaDrag = 0.02;
yDrag = 25;
yInertiaDrag = 8;
yWaveBreaking1 = 2; 
yWaveBreaking2 = 1.25;
x1 = xlim;
graphScale_x = x1(1);
legend show
legend('AutoUpdate','off')


 if x1(1) >= 0.01
     graphScale_x = 0.01;
       
 end

hold on
% deep water breaking limit (H/lambda_w = 0.14)
x2 = logspace(1,log10(graphScale_x),2);
y_breaking = 0.14*pi./x2;
plot(x2, y_breaking, 'k')
text(1,yWaveBreaking1,'{\it Wave Breaking}','HorizontalAlignment','center','FontSize',8)
text(1,yWaveBreaking2,'$\frac{H}{\lambda_w} > 0.14$','Interpreter','latex','HorizontalAlignment','center','FontSize',8)

hold on
gs_x = min(x2);
gs2 = max(x2);
% upper bound of low drag region
    ldv = 20;
    ldh = 0.14*pi/ldv;
    line([gs_x,ldh],[ldv,ldv], 'Color','k','LineStyle','--')
    text(xDrag,yDrag,'{\it Drag}','HorizontalAlignment','center','FontSize',8)

hold on
% upper bound of small drag region
    sdv = 1.5;
    sdh = 0.14*pi/sdv;
    line([gs_x,sdh],[sdv,sdv], 'Color','k','LineStyle','--')
    text(gs_x*2.5,yInertiaDrag,'{\it Inertia & Drag}','HorizontalAlignment','center','FontSize',8)

hold on        
% upper bound of negligible drag region
    ndv = 0.25;
    ndh = 0.14*pi/ndv;
    line([gs_x,ndh],[ndv,ndv], 'Color','k','LineStyle','--')
    text(8*gs_x,0.7,'{\it Large Inertia}','HorizontalAlignment','center','FontSize',8)

    text(8*gs_x,8e-2,'{\it All Inertia}','HorizontalAlignment','center','FontSize',8)
    
hold on
% left bound of diffraction region
    drh = 0.5;
    drv = 0.14*pi/drh;
    line([drh,drh],[gs_x,drv], 'Color','k','LineStyle','--')
    text(2,8e-2,'{\it Diffraction}','HorizontalAlignment','center','FontSize',8)

axis([gs_x gs2 min(y_breaking) max(y_breaking)])

len = strlength(options.savepath);
if len > 1
    saveas(figure, options.savepath);
end 

hold off
