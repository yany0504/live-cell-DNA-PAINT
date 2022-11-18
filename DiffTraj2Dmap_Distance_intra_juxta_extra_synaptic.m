%plot 2D heatmap (Diff coeff vs Traj) and Distance
%plot from Diff_Traj_Distance file.

%% Step 1: loading files and concatenate files together.
clc;clear;close all;
FilePath = 'E:\Data\DNA-PAINT(Live-cell)\SPT_tracking\cLTP_tracking\210225_livecell_DNA-PAINT_cLTP_GFP_or_antiHomer1-EYFP_DIV15\Coverslip3(antiHomer1-EYFP_uPAINT\Cell1\Results\';
OutputFile='20min';
levelStep=5; % 2D contour levelstep
intra=0.5;
juxta=2;
levelStep_syn=1;
levelStep_juxta=1;
x0 = [2,-1.5,1,-0.5,1,0.5,   1,-3,0.5,-1,0.5,0.2]; %Inital guess parameters for fitting
x0_diff = [x0(1),x0(2),x0(3),x0(7),x0(8),x0(9)];
x0_traj = [x0(1),x0(4),x0(5),x0(7),x0(10),x0(11)];
%% load
if ~exist('FileName1','var')|| isempty(FileName1)
    [userfilein, userdirin]=uigetfile({
         '*.xlsx','Data file (*.xlsx)';...
        '*.*','All Files (*.*)'},'Select the diffusion vs trajectory file to process',...
        FilePath, 'MultiSelect','on');
    FileName1=fullfile(userdirin,userfilein);
else
    if ~exist(FileName1,'file')
        fprintf('File not found: %s\n',FileName1);
        return;
    end
end

if iscell(FileName1);
    for i=1:length(FileName1);
        data1{i}=xlsread(FileName1{i});
    end
else
    data1=xlsread(FileName1);
end

if iscell(data1);
    DiffTrajDist=data1{1};
    for i=2:length(FileName1);
        DiffTrajDist=cat(1,DiffTrajDist,data1{i});
    end
else
    DiffTrajDist=data1;
end

%% Step 2: plot parameters
Diffbin=0.2;
Diffmin=-5.0;
Diffmax=0;

Trajbin=0.1;
Trajmin=-1.5;
Trajmax=0.5;

binsizeDist=0.05;
Distmin=0;
Distmax=2;

pos2D=[0.12,0.14,0.6,0.6];
posDiff=[0.12,0.75,0.6,0.2];
posTraj=[0.73, 0.14,0.2,0.6];

FigSize=[100 100 1000 800];
%% Step 3: Filter based on the range
Diffusion=DiffTrajDist(DiffTrajDist(:,4)>Diffmin,2);
Trajectory=DiffTrajDist(DiffTrajDist(:,4)>Diffmin,4);
Distance=DiffTrajDist(DiffTrajDist(:,4)>Diffmin,5);
ID=DiffTrajDist(DiffTrajDist(:,4)>Diffmin,6);
xlswrite([FilePath OutputFile '_selected_Diff_Traj_Dist.xlsx'], [Diffusion Trajectory Distance ID]);
%xlswrite([FilePath OutputFile '_Trajectory.xlsx'],Trajectory);
%xlswrite([FilePath OutputFile '_Distance.xlsx'], Distance);

%% Step 4: Distance measurement
figure;
Distance_filtered=Distance(Distance <=juxta);
hist_Distance=histogram(Distance_filtered,'BinWidth',0.01,'BinLimits',[Distmin,Distmax],'Normalization','probability');
ax=gca;
ax.YAxis.FontSize=25;
ax.XAxis.FontSize=25;
set(gca,'FontWeight','bold');
xlabel('Distance (\mum)','fontweight','bold','FontSize',25);
ylabel('Frequency','fontweight','bold','FontSize',25);
set(gcf, 'Position', FigSize);
xlim([0 2]);
saveas(gcf,[FilePath OutputFile '_Distance_pdf'],'fig');

figure;
hist_Distance=histogram(Distance_filtered,'BinWidth',0.01,'BinLimits',[Distmin,Distmax],'Normalization','cdf');
ax=gca;
ax.YAxis.FontSize=25;
ax.XAxis.FontSize=25;
set(gca,'FontWeight','bold');
xlabel('Distance (\mum)','fontweight','bold','FontSize',25);
ylabel('Cumulative Frequency','fontweight','bold','FontSize',25);
set(gcf, 'Position', FigSize);
ylim([0 1]);
xlim([0 2]);
saveas(gcf,[FilePath OutputFile '_Distance_cdf'],'fig');
%% Step 5: 2D map with contour

[N,c]=hist3([Diffusion Trajectory],'ctrs',{Diffmin:Diffbin:Diffmax Trajmin:Trajbin:Trajmax});
% [x,resnorm,residual,exitflag]=Fit2Dcontour(N,c,x0,0);

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
title('Total','fontweight','bold','FontSize',25);
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
saveas(gcf,[FilePath OutputFile '_2Dcontour'],'fig');





%% Step 6. distinguishing synaptic/juxtasynaptic/extrasynaptic
synaptic=DiffTrajDist(DiffTrajDist(:,5)<intra,:);
juxtasynaptic=DiffTrajDist(DiffTrajDist(:,5)>intra&DiffTrajDist(:,5)<juxta,:);
extrasynaptic=DiffTrajDist(DiffTrajDist(:,5)>juxta,:);

Diff_syn=synaptic(synaptic(:,4)>Diffmin,2);
Traj_syn=synaptic(synaptic(:,4)>Diffmin,4);
Dist_syn=synaptic(synaptic(:,4)>Diffmin,5);
ID_syn=synaptic(synaptic(:,4)>Diffmin,6);
xlswrite([FilePath OutputFile '_synaptic_Diff_Traj_Dist.xlsx'], [Diff_syn Traj_syn Dist_syn ID_syn]);

Diff_juxta=juxtasynaptic(juxtasynaptic(:,4)>Diffmin,2);
Traj_juxta=juxtasynaptic(juxtasynaptic(:,4)>Diffmin,4);
Dist_juxta=juxtasynaptic(juxtasynaptic(:,4)>Diffmin,5);
ID_juxta=juxtasynaptic(juxtasynaptic(:,4)>Diffmin,6);
xlswrite([FilePath OutputFile '_juxtasynaptic_Diff_Traj_Dist.xlsx'], [Diff_juxta Traj_juxta Dist_juxta ID_juxta]);

Diff_extra=extrasynaptic(extrasynaptic(:,4)>Diffmin,2);
Traj_extra=extrasynaptic(extrasynaptic(:,4)>Diffmin,4);
Dist_extra=extrasynaptic(extrasynaptic(:,4)>Diffmin,5);
ID_extra=extrasynaptic(extrasynaptic(:,4)>Diffmin,6);
xlswrite([FilePath OutputFile '_extrasynaptic_Diff_Traj_Dist.xlsx'], [Diff_extra Traj_extra Dist_extra ID_extra]);

%% 2D plot for intra- juxta- extra-synaptic.
% 
[N_syn,c_syn]=hist3([Diff_syn Traj_syn],'ctrs',{Diffmin:Diffbin:Diffmax Trajmin:Trajbin:Trajmax});
% [x_syn,resnorm_syn,residual_syn,exitflag_syn]=Fit2Dcontour(N_syn,c_syn,x0,0); % fitting with 2 Gaussian

figure;
subplot('Position',pos2D);
[M_syn,C_syn]=contourf(c_syn{1},c_syn{2},N_syn','-','LevelStep',levelStep_syn);
colormap jet
ax=gca;
ax.YAxis.FontSize=25;
ax.XAxis.FontSize=25;
set(gca,'FontWeight','bold');
xlabel('Diffusion coefficient (\mum^2/s, log)','fontweight','bold','FontSize',25);
ylabel('Trajectory range (\mum, log)','fontweight','bold','FontSize',25);

subplot('Position',posDiff);
hist_diff_syn=histogram(Diff_syn,'BinWidth',Diffbin,'BinLimits',[Diffmin-0.5*Diffbin,Diffmax+0.5*Diffbin],'Normalization','count');
% fit 2 Gaussian diffusion
[x_diff_syn,resnorm_diff_syn,residual_diff_syn,exitflag_diff_syn]=Fit1DHistograms(hist_diff_syn.Values,grid_diff,x0_diff);
fit_diff_syn=Two1DGaussFunction(x_diff_syn,grid_diff_hr);
xlim([Diffmin-0.5*Diffbin Diffmax+0.5*Diffbin])
hold on;
plot(grid_diff_hr,fit_diff_syn,'LineWidth',3);
hold off;

xlim([Diffmin-0.5*Diffbin Diffmax+0.5*Diffbin])
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Synaptic','fontweight','bold','FontSize',25);
subplot('Position',posTraj);
hist_traj_syn=histogram(Traj_syn,'BinWidth',Trajbin,'BinLimits',[Trajmin-0.5*Trajbin,Trajmax+0.5*Trajbin],'Orientation','horizontal','Normalization','count');
% fit 2 Gaussian Trajectory
[x_traj_syn,resnorm_traj_syn,residual_traj_syn,exitflag_traj_syn]=Fit1DHistograms(hist_traj_syn.Values,grid_traj,x0_traj);
fit_traj_syn=Two1DGaussFunction(x_traj_syn,grid_traj_hr);
hold on;
plot(fit_traj_syn,grid_traj_hr,'LineWidth',3);
hold off;

ylim([Trajmin-0.5*Trajbin Trajmax+0.5*Trajbin])
set(gca,'ytick',[]);
set(gca,'xtick',[]);
set(gcf, 'Position', FigSize);
saveas(gcf,[FilePath OutputFile '_synaptic_2Dcontour'],'fig');

% juxta synaptic
[N_juxta,c_juxta]=hist3([Diff_juxta Traj_juxta],'ctrs',{Diffmin:Diffbin:Diffmax Trajmin:Trajbin:Trajmax});
% [x_juxta,resnorm_juxta,residual_juxta,exitflag_juxta]=Fit2Dcontour(N_juxta,c_juxta,x0,0);

figure;
subplot('Position',pos2D);
[M_juxta,C_juxta]=contourf(c_juxta{1},c_juxta{2},N_juxta','-','LevelStep',levelStep_juxta);
colormap jet
ax=gca;
ax.YAxis.FontSize=25;
ax.XAxis.FontSize=25;
set(gca,'FontWeight','bold');
xlabel('Diffusion coefficient (\mum^2/s, log)','fontweight','bold','FontSize',25);
ylabel('Trajectory range (\mum, log)','fontweight','bold','FontSize',25);

subplot('Position',posDiff);
hist_diff_juxta=histogram(Diff_juxta,'BinWidth',Diffbin,'BinLimits',[Diffmin-0.5*Diffbin,Diffmax+0.5*Diffbin],'Normalization','count');

% fit 2 Gaussian diffusion
[x_diff_juxta,resnorm_diff_juxta,residual_diff_juxta,exitflag_diff_juxta]=Fit1DHistograms(hist_diff_juxta.Values,grid_diff,x0_diff);
fit_diff_juxta=Two1DGaussFunction(x_diff_juxta,grid_diff_hr);
xlim([Diffmin-0.5*Diffbin Diffmax+0.5*Diffbin])
hold on;
plot(grid_diff_hr,fit_diff_juxta,'LineWidth',3);
hold off;

xlim([Diffmin-0.5*Diffbin Diffmax+0.5*Diffbin])
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Juxta-synaptic','fontweight','bold','FontSize',25);
subplot('Position',posTraj);
hist_traj_juxta=histogram(Traj_juxta,'BinWidth',Trajbin,'BinLimits',[Trajmin-0.5*Trajbin,Trajmax+0.5*Trajbin],'Orientation','horizontal','Normalization','count');
% fit 2 Gaussian Trajectory
[x_traj_juxta,resnorm_traj_juxta,residual_traj_juxta,exitflag_traj_juxta]=Fit1DHistograms(hist_traj_juxta.Values,grid_traj,x0_traj);
fit_traj_juxta=Two1DGaussFunction(x_traj_juxta,grid_traj_hr);
hold on;
plot(fit_traj_juxta,grid_traj_hr,'LineWidth',3);
hold off;
ylim([Trajmin-0.5*Trajbin Trajmax+0.5*Trajbin])
set(gca,'ytick',[]);
set(gca,'xtick',[]);
set(gcf, 'Position', FigSize);
saveas(gcf,[FilePath OutputFile '_juxtasynaptic_2Dcontour'],'fig');


