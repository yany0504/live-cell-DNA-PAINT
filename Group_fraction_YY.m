function [Fraction,GMModel_thisdata]=Group_fraction_YY(Diffusion, Trajectory, GMModel,options);
% Diffmin,Diffbin,Diffmax,Trajmin,Trajmax,Trajbin, levelstep
arguments
    Diffusion double
    Trajectory double
    GMModel gmdistribution
    options.Diffrange double = [-4, -0.5, 0.2] % min, max, bin
    options.Trajrange double = [-1, 0.3, 0.1] % min, max, bin
    options.PlotStyle string = 'scatter'
    options.levelstep double = 10
    options.colorcode string = 'parula'
    options.LineColor string = 'none'
    options.title string = 'Title'
    options.text string = 'on'
    options.contourline string = '--w'
    options.contourlinewidth double = 1
    options.contourlevelstep double = 0.2
end
size_title=35;
size_label=28;



pos2D=[0.1,0.1,0.9,0.9];
FigSize=[100 100 1000 800];
colorcode=options.colorcode;
LineColor=options.LineColor;
% parameters
levelStep=options.levelstep;
Diffmin=options.Diffrange(1);
Diffmax=options.Diffrange(2);
Diffbin=options.Diffrange(3);
Trajmin=options.Trajrange(1);
Trajmax=options.Trajrange(2);
Trajbin=options.Trajrange(3);
plot_title=options.title;

opts=statset('MaxIter',10000);% max number iterations

grp=cluster(GMModel,[Diffusion,Trajectory]);
frac_grp1=sum((grp(:,1)==1))/length(grp);
frac_grp2=sum((grp(:,1)==2))/length(grp);
Fraction=[frac_grp1;frac_grp2];
gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(GMModel,[x0 y0]),x,y);
S=struct('mu',GMModel.mu,'Sigma',GMModel.Sigma,'ComponentProportion',[frac_grp1,frac_grp2]);
GMModel_thisdata=gmdistribution(GMModel.mu,GMModel.Sigma,Fraction);
gmPDF_thisdata = @(x,y) arrayfun(@(x0,y0) pdf(GMModel_thisdata,[x0 y0]),x,y);
figure;
switch options.PlotStyle
    case 'scatter'
        switch colorcode;
            case 'parula'
                colormap parula
                set(gca,'Color',[0.24 0.15 0.6]);
            case 'jet'
                colormap jet
                set(gca,'Color',[0 0 0.5]);
            case 'hsv'
                colormap hsv
                set(gca,'Color',[1 0 0]);
            case 'hot'
                colormap hot
                set(gca,'Color',[0 0 0]);
            case 'cool'
                colormap cool
                set(gca,'Color',[0 1 1]);
            case 'spring'
                colormap spring
                set(gca,'Color',[1 0 1]);
            case 'summer'
                colormap summer
                set(gca,'Color',[0 0.5 0.4]);
            case 'pink'
                colormap pink
                set(gca,'Color',[0.1 0 0]);
            case 'bone'
                colormap bone
                set(gca,'Color',[0 0 0]);
            case 'gray'
                colormap gray
                set(gca,'Color',[0 0 0]);
        end
        scatplot(Diffusion,Trajectory,'ms',10);
        colorbar('off')
        hold on;
        xlim([Diffmin Diffmax]);
        ylim([Trajmin Trajmax]);
        g=gca;
        fcontour(gmPDF_thisdata,[g.XLim, g.YLim],options.contourline,'LevelStep',options.contourlevelstep,'LineWidth',options.contourlinewidth);
        if options.text == 'on'
            text((GMModel.mu(1,1)-0.5),(GMModel.mu(1,2)+0.4),[num2str(100*round(frac_grp1,3)) '%' ],'Color','w','fontweight','bold','FontSize',size_label);
            text((GMModel.mu(2,1)-0.5),(GMModel.mu(2,2)+0.4),[num2str(100*round(frac_grp2,3)) '%' ],'Color','w','fontweight','bold','FontSize',size_label);
        else
        end
end
ax=gca;
ax.YAxis.FontSize=size_label;
ax.XAxis.FontSize=size_label;
set(gca,'FontWeight','bold');
xlabel('Diffusion coefficient (\mum^2/s, log)','fontweight','bold','FontSize',size_title);
ylabel('Trajectory range (\mum, log)','fontweight','bold','FontSize',size_title);
title(plot_title,'fontweight','bold','FontSize',size_title);
set(gcf, 'Position', FigSize);


end

