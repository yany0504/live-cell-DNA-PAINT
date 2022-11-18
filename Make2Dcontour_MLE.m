function [GMModel,GMModel_diff,GMModel_traj]=Make2Dcontour_MLE(Diffusion, Trajectory, options);
% Diffmin,Diffbin,Diffmax,Trajmin,Trajmax,Trajbin, levelstep
arguments
    Diffusion double
    Trajectory double
    options.Diffrange double = [-5, -0.5, 0.2] % min, max, bin
    options.Trajrange double = [-1.7, 0.3, 0.2] % min, max, bin
    options.PlotStyle string = 'scatter'
    options.levelstep double = 10
    options.colorcode string = 'parula'
    options.LineColor string = 'none'
    options.title string = 'Title'
    options.text string = 'on'
    options.histograms string = 'on'
    options.contour string = 'on'
    options.contourline string = '--w'
    options.contourlinewidth double = 1
    options.contourlevelstep double = 0.1
    options.Numpopulation double = 2
end
size_title=30;
size_label=25;
FigSize=[100 100 900 700];

if options.histograms == 'on'
    
    % figure setting
    pos2D=[0.12,0.14,0.6,0.6];
    posDiff=[0.12,0.75,0.6,0.2];
    posTraj=[0.73, 0.14,0.2,0.6];
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
    
    [N,c]=hist3([Diffusion Trajectory],'ctrs',{Diffmin:Diffbin:Diffmax Trajmin:Trajbin:Trajmax});
    GMModel=fitgmdist([Diffusion Trajectory],options.Numpopulation,'Options',opts);
    gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(GMModel,[x0 y0]),x,y);
    GMModel_diff=fitgmdist(Diffusion,options.Numpopulation,'Options',opts);
    gmPDF_diff=@(x) arrayfun(@(x0) pdf(GMModel_diff,x0),x);
    GMModel_traj=fitgmdist(Trajectory,options.Numpopulation,'Options',opts);
    gmPDF_traj=@(x) arrayfun(@(x0) pdf(GMModel_traj,x0),x);
    
    figure;
    subplot('Position',pos2D);
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
            if options.contour == 'on'
                fcontour(gmPDF,[g.XLim, g.YLim],options.contourline,'LevelStep',options.contourlevelstep,'LineWidth',options.contourlinewidth);
            else                
            end
            
            if options.text == 'on'
                for i=1:options.Numpopulation;
                    text((GMModel.mu(i,1)-0.5),(GMModel.mu(i,2)+0.4),[num2str(100*round(GMModel.ComponentProportion(i),3)) '%' ],'Color','w','fontweight','bold','FontSize',size_label);
                end
                %                 text((GMModel.mu(1,1)-0.5),(GMModel.mu(1,2)+0.4),[num2str(100*round(GMModel.ComponentProportion(1),3)) '%' ],'Color','w','fontweight','bold','FontSize',size_label);
                %                 text((GMModel.mu(2,1)-0.5),(GMModel.mu(2,2)+0.4),[num2str(100*round(GMModel.ComponentProportion(2),3)) '%' ],'Color','w','fontweight','bold','FontSize',size_label);
            else
            end
        case '2Dcontour'
            [M,C]=contourf(c{1},c{2},N','-','LevelStep',levelStep,'linecolor',LineColor);
            switch colorcode;
                case 'parula'
                    colormap parula
                case 'jet'
                    colormap jet
                case 'hsv'
                    colormap hsv
                case 'hot'
                    colormap hot
                case 'cool'
                    colormap cool
                case 'spring'
                    colormap spring
                case 'summer'
                    colormap summer
                case 'pink'
                    colormap pink
                case 'bone'
                    colormap bone
                case 'gray'
                    colormap gray
            end
            hold on;
            g=gca;
            fcontour(gmPDF,[g.XLim, g.YLim],options.contourline,'LevelStep',options.contourlevelstep,'LineWidth',options.contourlinewidth);
            if options.text == 'on'
                for i=1:options.Numpopulation;
                    text((GMModel.mu(i,1)-0.5),(GMModel.mu(i,2)+0.4),[num2str(100*round(GMModel.ComponentProportion(i),3)) '%' ],'Color','w','fontweight','bold','FontSize',size_label);
                end
                %                 text((GMModel.mu(1,1)-0.5),(GMModel.mu(1,2)+0.4),[num2str(100*round(GMModel.ComponentProportion(1),3)) '%' ],'Color','w','fontweight','bold','FontSize',size_label);
                %                 text((GMModel.mu(2,1)-0.5),(GMModel.mu(2,2)+0.4),[num2str(100*round(GMModel.ComponentProportion(2),3)) '%' ],'Color','w','fontweight','bold','FontSize',size_label);
            else
            end
            
            
    end
    
    hold off;
    ax=gca;
    ax.YAxis.FontSize=size_label;
    ax.XAxis.FontSize=size_label;
    set(gca,'FontWeight','bold');
    xlabel('Diffusion coefficient (\mum^2/s, log)','fontweight','bold','FontSize',size_title);
    ylabel('Trajectory range (\mum, log)','fontweight','bold','FontSize',size_title);
    
    % subplot for diffusion
    subplot('Position',posDiff);
    hist_diff=histogram(Diffusion,'BinWidth',Diffbin,'BinLimits',[Diffmin-0.5*Diffbin,Diffmax+0.5*Diffbin],'Normalization','pdf');
    xlim([Diffmin-0.5*Diffbin Diffmax+0.5*Diffbin])
    hold on;
    g2=gca;
    fplot(gmPDF_diff,[g2.XLim],'LineWidth',3);
    
    % lines for each gaussian function
    x_diff=linspace(Diffmin,Diffmax,300);
    gm_diff=cell(options.Numpopulation,1);
    gmPDF_diff=cell(options.Numpopulation,1);
    y_diff=cell(options.Numpopulation,1);
    for i=1:options.Numpopulation;
        gm_diff{i}=gmdistribution(GMModel_diff.mu(i),GMModel_diff.Sigma(i));
        gmPDF_diff{i}=@(x) arrayfun(@(x0) pdf(gm_diff{i},x0),x);
        y_diff{i}=GMModel_diff.ComponentProportion(1)*gmPDF_diff{i}(x_diff);
        plot(x_diff,y_diff{i},'LineStyle','--','LineWidth',2);
    end
    alpha(0.2)
    %
    %     gm_diff1 = gmdistribution(GMModel_diff.mu(1),GMModel_diff.Sigma(1));
    %     gmPDF_diff1=@(x) arrayfun(@(x0) pdf(gm_diff1,x0),x);
    %     y_diff1=GMModel_diff.ComponentProportion(1)*gmPDF_diff1(x_diff);
    %     gm_diff2 = gmdistribution(GMModel_diff.mu(2),GMModel_diff.Sigma(2));
    %     gmPDF_diff2=@(x) arrayfun(@(x0) pdf(gm_diff2,x0),x);
    %     y_diff2=GMModel_diff.ComponentProportion(2)*gmPDF_diff2(x_diff);
    %     plot(x_diff,y_diff1,'Color',[0.4660 0.6740 0.1880],'LineStyle','--','LineWidth',2);
    %     plot(x_diff,y_diff2,'Color',[0.3010 0.7450 0.9330],'LineStyle','--','LineWidth',2);
    %     alpha(0.2)
    
    if options.text == 'on'
        for i=1:options.Numpopulation;
            text((GMModel_diff.mu(i)-0.2),(max(hist_diff.Values)+0.1),[num2str(100*round(GMModel_diff.ComponentProportion(i),3)) '%' ],'Color','k','fontweight','bold','FontSize',size_label);
        end
%         text((GMModel_diff.mu(1)-0.2),(max(hist_diff.Values)+0.1),[num2str(100*round(GMModel_diff.ComponentProportion(1),3)) '%' ],'Color','k','fontweight','bold','FontSize',size_label);
%         text((GMModel_diff.mu(2)-0.2),(max(hist_diff.Values)+0.1),[num2str(100*round(GMModel_diff.ComponentProportion(2),3)) '%' ],'Color','k','fontweight','bold','FontSize',size_label);
    else
    end
    ylim([0 max(hist_diff.Values)+0.2]);
    hold off;
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    title(plot_title,'fontweight','bold','FontSize',size_title);
    % subplot for trajectory
    subplot('Position',posTraj);
    f=fplot(gmPDF_traj,[-2 2]);
    X=f.XData;
    Y=f.YData;
    hist_traj=histogram(Trajectory,'BinWidth',Trajbin,'BinLimits',[Trajmin-0.5*Trajbin,Trajmax+0.5*Trajbin],'Orientation','horizontal','Normalization','pdf');
    hold on;
    g3=gca;
    plot(Y,X,'LineWidth',3);
    
    % lines for each gaussian function
    x_traj=linspace(Trajmin,Trajmax,300);
    gm_traj=cell(options.Numpopulation,1);
    gmPDF_traj=cell(options.Numpopulation,1);
    y_traj=cell(options.Numpopulation,1);
    for i=1:options.Numpopulation;
        gm_traj{i}=gmdistribution(GMModel_traj.mu(i),GMModel_traj.Sigma(i));
        gmPDF_traj{i}=@(x) arrayfun(@(x0) pdf(gm_traj{i},x0),x);
        y_traj{i}=GMModel_traj.ComponentProportion(i)*gmPDF_traj{i}(x_traj);
        plot(y_traj{i},x_traj,'LineStyle','--','LineWidth',2);
    end
    alpha(0.2)
    %     gm_traj1 = gmdistribution(GMModel_traj.mu(1),GMModel_traj.Sigma(1));
    %     gmPDF_traj1=@(x) arrayfun(@(x0) pdf(gm_traj1,x0),x);
    %     y_traj1=GMModel_traj.ComponentProportion(1)*gmPDF_traj1(x_traj);
    %     gm_traj2 = gmdistribution(GMModel_traj.mu(2),GMModel_traj.Sigma(2));
    %     gmPDF_traj2=@(x) arrayfun(@(x0) pdf(gm_traj2,x0),x);
    %     y_traj2=GMModel_traj.ComponentProportion(2)*gmPDF_traj2(x_traj);
    %     plot(y_traj1,x_traj,'Color',[0.4660 0.6740 0.1880],'LineStyle','--','LineWidth',2);
    %     plot(y_traj2,x_traj,'Color',[0.3010 0.7450 0.9330],'LineStyle','--','LineWidth',2);
    %     alpha(0.2)
    
    if options.text == 'on'
        for i=1:options.Numpopulation;
            text((max(hist_traj.Values)+0.05),(GMModel_traj.mu(i)+0.1),[num2str(100*round(GMModel_traj.ComponentProportion(i),3)) '%' ],'Color','k','fontweight','bold','FontSize',size_label);
        end
%         text((max(hist_traj.Values)+0.05),(GMModel_traj.mu(1)+0.1),[num2str(100*round(GMModel_traj.ComponentProportion(1),3)) '%' ],'Color','k','fontweight','bold','FontSize',size_label);
%         text((max(hist_traj.Values)+0.05),(GMModel_traj.mu(2)+0.1),[num2str(100*round(GMModel_traj.ComponentProportion(2),3)) '%' ],'Color','k','fontweight','bold','FontSize',size_label);
    else
    end
    xlim([0 max(hist_traj.Values)+0.5]);
    hold off;
    ylim([Trajmin-0.5*Trajbin Trajmax+0.5*Trajbin])
    set(gca,'ytick',[]);
    set(gca,'xtick',[]);
    set(gcf, 'Position', FigSize);
    
    
elseif options.histograms == 'off'
    pos2D=[0.1,0.1,0.9,0.9];
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
    
    [N,c]=hist3([Diffusion Trajectory],'ctrs',{Diffmin:Diffbin:Diffmax Trajmin:Trajbin:Trajmax});
    GMModel=fitgmdist([Diffusion Trajectory],options.Numpopulation,'Options',opts);
    gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(GMModel,[x0 y0]),x,y);
    GMModel_diff=fitgmdist(Diffusion,options.Numpopulation,'Options',opts);
    gmPDF_diff=@(x) arrayfun(@(x0) pdf(GMModel_diff,x0),x);
    GMModel_traj=fitgmdist(Trajectory,options.Numpopulation,'Options',opts);
    gmPDF_traj=@(x) arrayfun(@(x0) pdf(GMModel_traj,x0),x);
    
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
            if options.contour == 'on'
                fcontour(gmPDF,[g.XLim, g.YLim],options.contourline,'LevelStep',options.contourlevelstep,'LineWidth',options.contourlinewidth);
            else
            end
            
            if options.text == 'on'
                for i=1:options.Numpopulation;
                    text((GMModel.mu(i,1)-0.5),(GMModel.mu(i,2)+0.4),[num2str(100*round(GMModel.ComponentProportion(i),3)) '%' ],'Color','w','fontweight','bold','FontSize',size_label);
                end
                %                 text((GMModel.mu(1,1)-0.5),(GMModel.mu(1,2)+0.4),[num2str(100*round(GMModel.ComponentProportion(1),3)) '%' ],'Color','w','fontweight','bold','FontSize',size_label);
                %                 text((GMModel.mu(2,1)-0.5),(GMModel.mu(2,2)+0.4),[num2str(100*round(GMModel.ComponentProportion(2),3)) '%' ],'Color','w','fontweight','bold','FontSize',size_label);
            else
            end
        case '2Dcontour'
            [M,C]=contourf(c{1},c{2},N','-','LevelStep',levelStep,'linecolor',LineColor);
            switch colorcode;
                case 'parula'
                    colormap parula
                case 'jet'
                    colormap jet
                case 'hsv'
                    colormap hsv
                case 'hot'
                    colormap hot
                case 'cool'
                    colormap cool
                case 'spring'
                    colormap spring
                case 'summer'
                    colormap summer
                case 'pink'
                    colormap pink
                case 'bone'
                    colormap bone
                case 'gray'
                    colormap gray
            end
            hold on;
            g=gca;
            fcontour(gmPDF,[g.XLim, g.YLim],options.contourline,'LevelStep',options.contourlevelstep,'LineWidth',options.contourlinewidth);
            if options.text == 'on'
                for i=1:options.Numpopulation;
                    text((GMModel.mu(i,1)-0.5),(GMModel.mu(i,2)+0.4),[num2str(100*round(GMModel.ComponentProportion(i),3)) '%' ],'Color','w','fontweight','bold','FontSize',size_label);
                end
                %                 text((GMModel.mu(1,1)-0.5),(GMModel.mu(1,2)+0.4),[num2str(100*round(GMModel.ComponentProportion(1),3)) '%' ],'Color','w','fontweight','bold','FontSize',size_label);
                %                 text((GMModel.mu(2,1)-0.5),(GMModel.mu(2,2)+0.4),[num2str(100*round(GMModel.ComponentProportion(2),3)) '%' ],'Color','w','fontweight','bold','FontSize',size_label);
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

end

