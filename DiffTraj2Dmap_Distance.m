%plot 2D heatmap (Diff coeff vs Traj) and Distance plot

%% Step 1: loading files and concatenate files together.
clc;clear;close all;
FilePath = 'E:\Data\DNA-PAINT(Live-cell)\SPT_tracking\cLTP_tracking\210225_livecell_DNA-PAINT_cLTP_GFP_or_antiHomer1-EYFP_DIV15\Coverslip3(antiHomer1-EYFP_uPAINT\Cell1\Results';
OutputFile='before';
levelStep=10; % 2D contour levelstep
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

if ~exist('FileName2','var')|| isempty(FileName2)
    [userfilein, userdirin]=uigetfile({
         '*.xlsx','Data file (*.xlsx)';...
        '*.*','All Files (*.*)'},'Select the distance file to process',...
        FilePath, 'MultiSelect','on');
    FileName2=fullfile(userdirin,userfilein);
else
    if ~exist(FileName2,'file')
        fprintf('File not found: %s\n',FileName2);
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
if iscell(FileName2);
    for i=1:length(FileName2);
        data2{i}=xlsread(FileName2{i});
    end
else
    data2=xlsread(FileName2);
end

if iscell(data1);
    DiffTraj=data1{1};
    for i=2:length(FileName1);
        DiffTraj=cat(1,DiffTraj,data1{i});
    end
else
    DiffTraj=data1;
end
if iscell(data2);
    Distance=data2{1};
    for i=2:length(FileName2);
        Distance=cat(1,Distance,data2{i});
    end
else
    Distance=data2;
end

%% Step 2: plot parameters
Diffbin=0.2;
Diffmin=-6.0;
Diffmax=0;

Trajbin=0.1;
Trajmin=-1.5;
Trajmax=0.5;

binsizeDist=0.05;
Distmin=0;
Distmax=2;

pos2D=[0.15,0.15,0.6,0.6];
posDiff=[0.15,0.77,0.6,0.2];
posTraj=[0.77, 0.15,0.2,0.6];

FigSize=[100 100 1000 800];
%% Step 3: Filter based on the range
Diffusion=DiffTraj(DiffTraj(:,4)>Diffmin,4);
Trajectory=DiffTraj(DiffTraj(:,4)>Diffmin,3);
xlswrite([FilePath OutputFile '_Diffusion.xlsx'], Diffusion);
xlswrite([FilePath OutputFile '_Trajectory.xlsx'],Trajectory);
xlswrite([FilePath OutputFile '_Distance.xlsx'], Distance);
%% Step 3.5: 2D map with contour
%Diffusion=DiffTraj(:,4);
%Trajectory=DiffTraj(:,3);
[N,c]=hist3([Diffusion Trajectory],'ctrs',{Diffmin:Diffbin:Diffmax Trajmin:Trajbin:Trajmax});
figure;
subplot('Position',pos2D);
[M,C]=contourf(c{1},c{2},N','-','LevelStep',levelStep);
colormap jet
ax=gca;
ax.YAxis.FontSize=25;
ax.XAxis.FontSize=25;
set(gca,'FontWeight','bold');
xlabel('Diffusion coefficient (\mum^2/s, log)','fontweight','bold','FontSize',25);
ylabel('Trajectory length (\mum, log)','fontweight','bold','FontSize',25);

subplot('Position',posDiff);
hist_diff=histogram(DiffTraj(:,4),'BinWidth',Diffbin,'BinLimits',[Diffmin-0.5*Diffbin,Diffmax+0.5*Diffbin],'Normalization','probability');
xlim([Diffmin-0.5*Diffbin Diffmax+0.5*Diffbin])
set(gca,'xtick',[]);
set(gca,'ytick',[]);
subplot('Position',posTraj);
hist_traj=histogram(DiffTraj(:,3),'BinWidth',Trajbin,'BinLimits',[Trajmin-0.5*Trajbin,Trajmax+0.5*Trajbin],'Orientation','horizontal','Normalization','probability');
ylim([Trajmin-0.5*Trajbin Trajmax+0.5*Trajbin])
set(gca,'ytick',[]);
set(gca,'xtick',[]);
set(gcf, 'Position', FigSize);
saveas(gcf,[FilePath OutputFile '_2Dcontour'],'fig');


%% Step 4: Distance measurement
figure;
hist_Distance=histogram(Distance,'BinWidth',0.01,'BinLimits',[Distmin,Distmax],'Normalization','probability');
ax=gca;
ax.YAxis.FontSize=25;
ax.XAxis.FontSize=25;
set(gca,'FontWeight','bold');
xlabel('Distance (\mum)','fontweight','bold','FontSize',25);
ylabel('Frequency','fontweight','bold','FontSize',25);
set(gcf, 'Position', FigSize);

saveas(gcf,[FilePath OutputFile '_Distance_pdf'],'fig');
figure;
hist_Distance=histogram(Distance,'BinWidth',0.01,'BinLimits',[Distmin,Distmax],'Normalization','cdf');
ax=gca;
ax.YAxis.FontSize=25;
ax.XAxis.FontSize=25;
set(gca,'FontWeight','bold');
xlabel('Distance (\mum)','fontweight','bold','FontSize',25);
ylabel('Cumulative Frequency','fontweight','bold','FontSize',25);
set(gcf, 'Position', FigSize);
ylim([0 1]);
saveas(gcf,[FilePath OutputFile '_Distance_cdf'],'fig');

