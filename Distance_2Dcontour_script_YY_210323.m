%plot 2D heatmap (Diff coeff vs Traj) and Distance
%plot from Diff_Traj_Distance file.

%% Step 1: loading files and concatenate files together.
clc;clear;close all;
FilePath = 'E:\Data\DNA-PAINT(Live-cell)\SPT_tracking\Basal_tracking\210304_DNA-PAINT_1xR3_antiHomer1-mGeos_3D_DIV15_HP\7xR3_1nM_R3_9nt-LD655\Results\';
OutputFile='Testout_before';
levelStep=10; % 2D contour levelstep
intra=0.5;
juxta=2;
levelStep_syn=5;
levelStep_juxta=5;
x0 = [2,-1.5,1,-0.5,1,0.5,   1,-3,0.5,-1,0.5,0.2]; %Inital guess parameters for fitting
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
[x,x_diff,x_traj]=Make2Dcontour(Diffusion,Trajectory,x0,'levelstep',levelStep,'title','Total');

%saveas(gcf,[FilePath OutputFile '_2Dcontour'],'fig');


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
[x_syn,x_diff_syn,x_traj_syn]=Make2Dcontour(Diff_syn,Traj_syn,x0,'levelstep',levelStep_syn,'title','Synaptic');
N_syn_diff_1=sqrt(2*pi)*x_diff_syn(1)*x_diff_syn(3);
N_syn_diff_2=sqrt(2*pi)*x_diff_syn(4)*x_diff_syn(6);
N_syn_traj_1=sqrt(2*pi)*x_traj_syn(1)*x_traj_syn(3);
N_syn_traj_2=sqrt(2*pi)*x_traj_syn(4)*x_traj_syn(6);
N_syn_1=2*pi*x_syn(1)*x_syn(3)*x_syn(5);
N_syn_2=2*pi*x_syn(7)*x_syn(9)*x_syn(11);
[x_juxta,x_diff_juxta,x_traj_juxta]=Make2Dcontour(Diff_juxta,Traj_juxta,x0,'levelstep',levelStep_juxta,'title','Juxta-synaptic');
N_juxta_diff_1=sqrt(2*pi)*x_diff_juxta(1)*x_diff_juxta(3);
N_juxta_diff_2=sqrt(2*pi)*x_diff_juxta(4)*x_diff_juxta(6);
N_juxta_traj_1=sqrt(2*pi)*x_traj_juxta(1)*x_traj_juxta(3);
N_juxta_traj_2=sqrt(2*pi)*x_traj_juxta(4)*x_traj_juxta(6);
N_juxta_1=2*pi*x_juxta(1)*x_juxta(3)*x_juxta(5);
N_juxta_2=2*pi*x_juxta(7)*x_juxta(9)*x_juxta(11);