function [x,x_diff,x_traj]=Make2Dcontour(Diffusion, Trajectory, x0, options);
% Diffmin,Diffbin,Diffmax,Trajmin,Trajmax,Trajbin, levelstep
arguments
    Diffusion double
    Trajectory double
    x0 double
    options.Diffrange double = [-5, 0, 0.2] % min, max, bin
    options.Trajrange double = [-1.5, 0.5, 0.1] % min, max, bin
    options.levelstep double = 10
    options.title string = 'Title'
    
end
% figure setting
pos2D=[0.12,0.14,0.6,0.6];
posDiff=[0.12,0.75,0.6,0.2];
posTraj=[0.73, 0.14,0.2,0.6];

FigSize=[100 100 1000 800];
% parameters
levelStep=options.levelstep;
Diffmin=options.Diffrange(1);
Diffmax=options.Diffrange(2);
Diffbin=options.Diffrange(3);
Trajmin=options.Trajrange(1);
Trajmax=options.Trajrange(2);
Trajbin=options.Trajrange(3);
x0_diff = [x0(1),x0(2),x0(3),x0(7),x0(8),x0(9)];
x0_traj = [x0(1),x0(4),x0(5),x0(7),x0(10),x0(11)];
plot_title=options.title;

[N,c]=hist3([Diffusion Trajectory],'ctrs',{Diffmin:Diffbin:Diffmax Trajmin:Trajbin:Trajmax});
[x,resnorm,residual,exitflag]=Fit2Dcontour(N,c,x0,0);

figure;
subplot('Position',pos2D);
[M,C]=contourf(c{1},c{2},N','-','LevelStep',levelStep);
colormap jet
ax=gca;
ax.YAxis.FontSize=25;
ax.XAxis.FontSize=25;
set(gca,'FontWeight','bold');
xlabel('Diffusion coefficient (\mum^2/s, log)','fontweight','bold','FontSize',25);
ylabel('Trajectory range (\mum, log)','fontweight','bold','FontSize',25);

subplot('Position',posDiff);
hist_diff=histogram(Diffusion,'BinWidth',Diffbin,'BinLimits',[Diffmin-0.5*Diffbin,Diffmax+0.5*Diffbin],'Normalization','count');
% fit 2 Gaussian Diffusion
grid_diff=linspace(Diffmin,Diffmax,hist_diff.NumBins);
grid_diff_hr=linspace(Diffmin,Diffmax,100);
[x_diff,resnorm_diff,residual_diff,exitflag_diff]=Fit1DHistograms(hist_diff.Values,grid_diff,x0_diff);
fit_diff=Two1DGaussFunction(x_diff,grid_diff_hr);
xlim([Diffmin-0.5*Diffbin Diffmax+0.5*Diffbin])
hold on;
plot(grid_diff_hr,fit_diff,'LineWidth',3);
hold off;

set(gca,'xtick',[]);
set(gca,'ytick',[]);
title(plot_title,'fontweight','bold','FontSize',25);
subplot('Position',posTraj);
hist_traj=histogram(Trajectory,'BinWidth',Trajbin,'BinLimits',[Trajmin-0.5*Trajbin,Trajmax+0.5*Trajbin],'Orientation','horizontal','Normalization','count');
% fit 2 Gaussian Trajectory
grid_traj=linspace(Trajmin,Trajmax,hist_traj.NumBins);
grid_traj_hr=linspace(Trajmin,Trajmax,100);
[x_traj,resnorm_traj,residual_traj,exitflag_traj]=Fit1DHistograms(hist_traj.Values,grid_traj,x0_traj);
fit_traj=Two1DGaussFunction(x_traj,grid_traj_hr);
hold on;
plot(fit_traj,grid_traj_hr,'LineWidth',3);
hold off;

ylim([Trajmin-0.5*Trajbin Trajmax+0.5*Trajbin])
set(gca,'ytick',[]);
set(gca,'xtick',[]);
set(gcf, 'Position', FigSize);

end
