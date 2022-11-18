%plot 2D heatmap (Diff coeff vs Traj)
%plot from Diff_Traj_Distance file.

%% Step 1: loading files and concatenate files together.
clc;clear;%close all;
FilePath = 'H:\Data_from07222021\DNA-PAINT(Live-cell)\SPT_tracking\cLTP_tracking\211202_cLTP_DNA-PAINT_LiefAct_DIV15\Results\';
OutputFile='after';
levelStep=.1; % 2D contour levelstep


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
Diffmin=-4.5;
Diffmax=-0.5;
Trajbin=0.1;
Trajmin=-1.5;
Trajmax=0.2;
binsizeDist=0.05;
Distmin=0;
Distmax=2;

pos2D=[0.12,0.14,0.6,0.6];
posDiff=[0.12,0.75,0.6,0.2];
posTraj=[0.73, 0.14,0.2,0.6];

FigSize=[100 100 1000 800];
%% Step 3: Filter based on the range
Diffusion=DiffTrajDist(:,4);
Trajectory=DiffTrajDist(:,3);

%% 2D plot for total.

[GMModel_tot,GMModel_diff_tot,GMModel_traj_tot]=Make2Dcontour_MLE(Diffusion,Trajectory,'Diffrange',[-4.5, -0.5, 0.2],'Trajrange',[-1.5, 0.3, 0.1],'contourlevelstep',levelStep,'title','Total',...
    'PlotStyle','scatter','LineColor','none','colorcode','parula','text','on','histograms','off','contour','on','Numpopulation',2);
saveas(gcf,[FilePath OutputFile '_total_2D'],'fig');

%% cluster analysis
grp_tot=cluster(GMModel_tot,[Diffusion,Trajectory]);

figure;
scatter(Diffusion,Trajectory,20,grp_tot);
xlim([-6 0]);ylim([-2 0.5]);
title('Total');
