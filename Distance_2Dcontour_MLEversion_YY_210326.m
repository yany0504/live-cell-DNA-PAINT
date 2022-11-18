%plot 2D heatmap (Diff coeff vs Traj) and Distance
%plot from Diff_Traj_Distance file.

%% Step 1: loading files and concatenate files together.
clc;clear;
close all;
FilePath = 'H:\Data_from07222021\DNA-PAINT(Live-cell)\SPT_tracking\Basal_tracking\211217_basalTracking_37C\Coverslip3_antiHomer1_5xR1_SPT3D\Results\';
OutputFile='all';
levelStep=1; % 2D contour levelstep
levelStep_syn=1;
levelStep_juxta=1;

intra=0.3;
juxta=2.0;
plot_Diffrange = [-4.5, -0.5, 0.2];
plot_Trajrange = [-1.5, 0.3, 0.2];

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
Diffmin=-4;
Diffmax=-0.5;
Trajbin=0.1;
Trajmin=-1.5;
Trajmax=0.2;
binsizeDist=0.02;
Distmin=0;
Distmax=2;

pos2D=[0.12,0.14,0.6,0.6];
posDiff=[0.12,0.75,0.6,0.2];
posTraj=[0.73, 0.14,0.2,0.6];

FigSize=[100 100 1200 900];
%% Step 3: Filter based on the range
Diffusion=DiffTrajDist(DiffTrajDist(:,4)>Diffmin,2);
Trajectory=DiffTrajDist(DiffTrajDist(:,4)>Diffmin,4);
Distance=DiffTrajDist(DiffTrajDist(:,4)>Diffmin,5);
ID=DiffTrajDist(DiffTrajDist(:,4)>Diffmin,6);
if size(DiffTrajDist,2) >6
    frame=DiffTrajDist(DiffTrajDist(:,4)>Diffmin,7);
    synapse=DiffTrajDist(DiffTrajDist(:,4)>Diffmin,8);
    xlswrite([FilePath OutputFile '_selected_Diff_Traj_Dist.xlsx'], [Diffusion Trajectory Distance ID frame synapse]);
else
    xlswrite([FilePath OutputFile '_selected_Diff_Traj_Dist.xlsx'], [Diffusion Trajectory Distance ID]);
end



%% Step 4: Distance measurement
figure;
Distance_filtered=Distance(Distance <=juxta);
hist_Distance=histogram(Distance_filtered,'BinWidth',binsizeDist,'BinLimits',[Distmin,Distmax],'Normalization','probability');
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
hist_Distance=histogram(Distance_filtered,'BinWidth',binsizeDist,'BinLimits',[Distmin,Distmax],'Normalization','cdf');
ax=gca;
ax.YAxis.FontSize=30;
ax.XAxis.FontSize=30;
set(gca,'FontWeight','bold');
xlabel('Distance (\mum)','fontweight','bold','FontSize',35);
ylabel('Cumulative Frequency','fontweight','bold','FontSize',35);
set(gcf, 'Position', FigSize);
ylim([0 1]);
xlim([0 2]);
saveas(gcf,[FilePath OutputFile '_Distance_cdf'],'fig');
%% Step 5. distinguishing synaptic/juxtasynaptic/extrasynaptic
synaptic=DiffTrajDist(DiffTrajDist(:,5)<intra,:);
juxtasynaptic=DiffTrajDist(DiffTrajDist(:,5)>intra&DiffTrajDist(:,5)<juxta,:);
extrasynaptic=DiffTrajDist(DiffTrajDist(:,5)>juxta,:);

Diff_syn=synaptic(synaptic(:,4)>Diffmin,2);
Traj_syn=synaptic(synaptic(:,4)>Diffmin,4);
Dist_syn=synaptic(synaptic(:,4)>Diffmin,5);
ID_syn=synaptic(synaptic(:,4)>Diffmin,6);
if size(DiffTrajDist,2) >6
    frame_syn=synaptic(synaptic(:,4)>Diffmin,7);
    synapse_syn=synaptic(synaptic(:,4)>Diffmin,8);
    xlswrite([FilePath OutputFile '_synaptic_Diff_Traj_Dist.xlsx'], [Diff_syn Traj_syn Dist_syn ID_syn frame_syn synapse_syn]);
else
    xlswrite([FilePath OutputFile '_synaptic_Diff_Traj_Dist.xlsx'], [Diff_syn Traj_syn Dist_syn ID_syn]);
end

Diff_juxta=juxtasynaptic(juxtasynaptic(:,4)>Diffmin,2);
Traj_juxta=juxtasynaptic(juxtasynaptic(:,4)>Diffmin,4);
Dist_juxta=juxtasynaptic(juxtasynaptic(:,4)>Diffmin,5);
ID_juxta=juxtasynaptic(juxtasynaptic(:,4)>Diffmin,6);
if size(DiffTrajDist,2)>6
    frame_juxta=juxtasynaptic(juxtasynaptic(:,4)>Diffmin,7);
    synapse_juxta=juxtasynaptic(juxtasynaptic(:,4)>Diffmin,8);
    xlswrite([FilePath OutputFile '_juxtasynaptic_Diff_Traj_Dist.xlsx'], [Diff_juxta Traj_juxta Dist_juxta ID_juxta frame_juxta synapse_juxta]);
else
    xlswrite([FilePath OutputFile '_juxtasynaptic_Diff_Traj_Dist.xlsx'], [Diff_juxta Traj_juxta Dist_juxta ID_juxta]);
end

Diff_extra=extrasynaptic(extrasynaptic(:,4)>Diffmin,2);
Traj_extra=extrasynaptic(extrasynaptic(:,4)>Diffmin,4);
Dist_extra=extrasynaptic(extrasynaptic(:,4)>Diffmin,5);
ID_extra=extrasynaptic(extrasynaptic(:,4)>Diffmin,6);
xlswrite([FilePath OutputFile '_extrasynaptic_Diff_Traj_Dist.xlsx'], [Diff_extra Traj_extra Dist_extra ID_extra]);

Diff_tot=[Diff_syn; Diff_juxta];
Traj_tot=[Traj_syn; Traj_juxta];
Dist_tot=[Dist_syn; Dist_juxta];
ID_tot=[ID_syn; ID_juxta];
%% 2D plot for synaptic, juxta-, total.
levelStep=20;
levelStep_syn=10;
levelStep_juxta=0.1;
[GMModel_tot,GMModel_diff_tot,GMModel_traj_tot]=Make2Dcontour_MLE(Diff_tot,Traj_tot,'levelstep',levelStep,'title','Total',...
    'PlotStyle','scatter','LineColor','none','colorcode','parula','text','on','histograms','off','contour','on','Numpopulation',2,'Diffrange',plot_Diffrange,'Trajrange',plot_Trajrange);
saveas(gcf,[FilePath OutputFile '_total_2D'],'fig');

% [Fraction_synaptic,GMModel_syn]=Group_fraction_YY(Diff_syn,Traj_syn,GMModel_tot,'title','Synaptic');
% saveas(gcf,[FilePath OutputFile '_syn_2D_usingGMModeltotal'],'fig');
% 
[GMModel_syn,GMModel_diff_syn,GMModel_traj_syn]=Make2Dcontour_MLE(Diff_syn,Traj_syn,'levelstep',levelStep_syn,'title','Synaptic',...
    'PlotStyle','scatter','LineColor','none','colorcode','parula','text','on','histograms','off','contour','on','Numpopulation',2,'Diffrange',plot_Diffrange,'Trajrange',plot_Trajrange);
saveas(gcf,[FilePath OutputFile '_syn_2D_usingitself'],'fig');

% [Fraction_juxta,GMModel_juxta]=Group_fraction_YY(Diff_juxta,Traj_juxta,GMModel_tot,'title','Juxta-synaptic');
% saveas(gcf,[FilePath OutputFile '_juxta_2D_usingGMModeltotal'],'fig');

[GMModel_juxt,GMModel_diff_juxt,GMModel_traj_juxta]=Make2Dcontour_MLE(Diff_juxta,Traj_juxta,'levelstep',levelStep_juxta,'title','Juxta-synaptic',...
   'PlotStyle','scatter','LineColor','none','colorcode','parula','text','on','histograms','off','contour','on','Numpopulation',2,'Diffrange',plot_Diffrange,'Trajrange',plot_Trajrange);
saveas(gcf,[FilePath OutputFile '_juxta_2D_usingitself'],'fig');

save([FilePath OutputFile '_GMModel_analysis.mat']);

%% fraction of immobile at different region.
Diff_dist=cell(20,1);
Traj_dist=cell(20,1);
grp_dist=cell(20,1);
count_immo_mobile=zeros(20,2);

for i=1:20
    Diff_dist{i}=DiffTrajDist(DiffTrajDist(:,5)<0.1*i&DiffTrajDist(:,5)>0.1*(i-1),2);
    Traj_dist{i}=DiffTrajDist(DiffTrajDist(:,5)<0.1*i&DiffTrajDist(:,5)>0.1*(i-1),4);
    grp_dist{i}=cluster(GMModel_tot,[Diff_dist{i},Traj_dist{i}]);
    count_immo_mobile(i,1)=sum(grp_dist{i}==1);
    count_immo_mobile(i,2)=sum(grp_dist{i}==2);
end

total_dist=sum(count_immo_mobile,2);
ratio_immobile=count_immo_mobile(:,1)./total_dist;

%% histogram for diffusion and trajectory
figure;
set(gcf,'Position',[100 100 1000 800]);
t=tiledlayout(2,2,'TileSpacing','Compact');

nexttile
h_diff_syn=histogram(Diff_syn,'BinWidth',0.2);
xlim([-5.5 -0.5])
set(gca,'FontSize',15)
ylabel('Count','FontSize',25)
title('Synaptic','FontSize',25)
ax = gca;
ax.TitleHorizontalAlignment = 'left';

nexttile
h_traj_syn=histogram(Traj_syn,'BinWidth',0.1);
xlim([-1.7 0.5])
set(gca,'FontSize',15)
nexttile
title('Juxta-synaptic')
h_diff_juxt=histogram(Diff_juxta,'BinWidth',0.2);
xlim([-5.5 -0.5])
set(gca,'FontSize',15)
xlabel('Diff. Coeff. (\mum^2/s, log)','FontSize',25)
ylabel('Count','FontSize',25)
title('Juxta-ynaptic','FontSize',25)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
nexttile
h_traj_juxta=histogram(Traj_juxta,'BinWidth',0.1);
xlim([-1.7 0.5])
set(gca,'FontSize',15)
xlabel('Traj range (\mum, log)','FontSize',25)

%% cluster analysis with Gaussian mixture model
% grp_syn=cluster(GMModel_tot,[Diff_syn,Traj_syn]);
% grp_juxta=cluster(GMModel_tot,[Diff_juxta,Traj_juxta]);
% % %grp_tot=cluster(GMModel_tot,[Diff_tot,Traj_tot]);
% % 
% figure;
% gscatter(Diff_syn,Traj_syn,grp_syn,'bgm');
% Diff_syn1=Diff_syn(grp_syn==1);
% Traj_syn1=Traj_syn(grp_syn==1);
% Diff_syn2=Diff_syn(grp_syn==2);
% Traj_syn2=Traj_syn(grp_syn==2);
% 
% Diff_juxta1=Diff_juxta(grp_juxta==1);
% Traj_juxta1=Traj_juxta(grp_juxta==1);
% Diff_juxta2=Diff_juxta(grp_juxta==2);
% Traj_juxta2=Traj_juxta(grp_juxta==2);

% ylabel('Trajectory range (log, \mum)');
% xlabel('Diffusion coefficient (log, /mum^2/s)');
% title('Synaptic GMM');
% xlim([-6 0]);ylim([-2 0.5]);
% 
% figure;
% gscatter(Diff_juxta,Traj_juxta,grp_juxta,'bgm');
% title('Juxtasynaptic');
% xlim([-6 0]);ylim([-2 0.5]);
% xlabel('Diffusion coefficient (log, /mum^2/s)');
% % 
% % subplot(1,3,3);scatter(Diff_tot,Traj_tot,20,grp_tot);
% % xlim([-6 0]);ylim([-2 0.5]);
% % title('Total');
% % 
% % set(gcf,'Position',[100 100 1400 400])

%% part 2
% making plot of diffusion coefficient as a function of distance from
% Homer1
% figure;
% Output=scatplot(Dist_tot,Diff_tot);
% xlabel('Distance from NN-Homer1 (\mum)');
% ylabel('Diffusion coefficient (log, \mum^2/s)');
% title('Total')
FigSize=[100, 100, 800, 700];
figure;
Output=scatplot(Diff_tot,Dist_tot);
ylabel('Distance from NN-Homer1 (\mum)','FontSize',30,'fontweight','Bold');
xlabel('Diffusion coefficient (\mum^2/s, log)','FontSize',30,'fontweight','Bold');
% title('Total')
set(gca,'Color',[0.24 0.15 0.6],'FontSize',25);
xlim([-5.2 -0.5]);
ylim([-0.05 2.05]);
syn_immo=Diff_tot(Diff_tot(:,1)<-2.5&Dist_tot<0.4);
set(gcf, 'Position', FigSize);
saveas(gcf,[FilePath OutputFile '_Diff_Distance_plot'],'fig');
% figure;
% Output=scatplot(Dist_syn,Diff_syn);
% xlabel('Distance from NN-Homer1 (\mum)');
% ylabel('Diffusion coefficient (log, \mum^2/s)');
% title('Synaptic')
% figure;
% Output=scatplot(Diff_syn,Dist_syn);
% ylabel('Distance from NN-Homer1 (\mum)');
% xlabel('Diffusion coefficient (log, \mum^2/s)');
% title('Synaptic')
% set(gca,'Color',[0.24 0.15 0.6]);
% xlim([-5.5 -0.5]);
% figure;
% Output=scatplot(Dist_juxta,Diff_juxta);
% xlabel('Distance from NN-Homer1 (\mum)');
% ylabel('Diffusion coefficient (log, \mum^2/s)');
% title('Juxta-synaptic')
% figure;
% Output=scatplot(Diff_juxta,Dist_juxta);
% ylabel('Distance from NN-Homer1 (\mum)');
% xlabel('Diffusion coefficient (log, \mum^2/s)');
% title('Juxta-synaptic')
% set(gca,'Color',[0.24 0.15 0.6]);
% xlim([-5.5 -0.5]);


dist_diff=[Dist_tot, Diff_tot];



diff_100=dist_diff(dist_diff(:,1)<=0.1,2);
diff_200=dist_diff(dist_diff(:,1)<=0.2&dist_diff(:,1)>0.1,2);
diff_300=dist_diff(dist_diff(:,1)<=0.3&dist_diff(:,1)>0.2,2);
diff_400=dist_diff(dist_diff(:,1)<=0.4&dist_diff(:,1)>0.3,2);
diff_500=dist_diff(dist_diff(:,1)<=0.5&dist_diff(:,1)>0.4,2);
diff_600=dist_diff(dist_diff(:,1)<=0.6&dist_diff(:,1)>0.5,2);
diff_700=dist_diff(dist_diff(:,1)<=0.7&dist_diff(:,1)>0.6,2);
diff_800=dist_diff(dist_diff(:,1)<=0.8&dist_diff(:,1)>0.7,2);
diff_900=dist_diff(dist_diff(:,1)<=0.9&dist_diff(:,1)>0.8,2);
diff_1000=dist_diff(dist_diff(:,1)<=1.0&dist_diff(:,1)>0.9,2);
diff_1100=dist_diff(dist_diff(:,1)<=1.1&dist_diff(:,1)>1.0,2);
diff_1200=dist_diff(dist_diff(:,1)<=1.2&dist_diff(:,1)>1.1,2);
diff_1300=dist_diff(dist_diff(:,1)<=1.3&dist_diff(:,1)>1.2,2);
diff_1400=dist_diff(dist_diff(:,1)<=1.4&dist_diff(:,1)>1.3,2);
diff_1500=dist_diff(dist_diff(:,1)<=1.5&dist_diff(:,1)>1.4,2);
diff_1600=dist_diff(dist_diff(:,1)<=1.6&dist_diff(:,1)>1.5,2);
diff_1700=dist_diff(dist_diff(:,1)<=1.7&dist_diff(:,1)>1.6,2);
diff_1800=dist_diff(dist_diff(:,1)<=1.8&dist_diff(:,1)>1.7,2);
diff_1900=dist_diff(dist_diff(:,1)<=1.9&dist_diff(:,1)>1.8,2);
diff_2000=dist_diff(dist_diff(:,1)<=2.0&dist_diff(:,1)>1.9,2);


%% calculate immobile/total
criteria=-2.5;
x_dist=[0.1:0.1:2];
y_immo_frac=[length(diff_100(diff_100<criteria))/length(diff_100),length(diff_200(diff_200<criteria))/length(diff_200),length(diff_300(diff_300<criteria))/length(diff_300),length(diff_400(diff_400<criteria))/length(diff_400),...
    length(diff_500(diff_500<criteria))/length(diff_500),length(diff_600(diff_600<criteria))/length(diff_600),length(diff_700(diff_700<criteria))/length(diff_700),length(diff_800(diff_800<criteria))/length(diff_800),...
    length(diff_900(diff_900<criteria))/length(diff_900),length(diff_1000(diff_1000<criteria))/length(diff_1000),length(diff_1100(diff_1100<criteria))/length(diff_1100),length(diff_1200(diff_1200<criteria))/length(diff_1200),...
    length(diff_1300(diff_1300<criteria))/length(diff_1300),length(diff_1400(diff_1400<criteria))/length(diff_1400),length(diff_1500(diff_1500<criteria))/length(diff_1500),length(diff_1600(diff_1600<criteria))/length(diff_1600),...
    length(diff_1700(diff_1700<criteria))/length(diff_1700),length(diff_1800(diff_1800<criteria))/length(diff_1800),length(diff_1900(diff_1900<criteria))/length(diff_1900),length(diff_2000(diff_2000<criteria))/length(diff_2000)];
figure;

plot(x_dist,y_immo_frac);
set(gca,'FontSize',20);
xlabel('Distance from NN-Homer1 (\mum)','FontSize',25);
ylabel('Fraction immobile GluA1-AMPAR','FontSize',25);
set(gcf, 'Position', FigSize);


%% plot every window
% 
% figure;
% h100=histogram(diff_100,'BinWidth',0.2,'Normalization','Probability');
% xlim([-6 -0.5]);
% xlabel('log(diffusion coefficient) (\mum^2/s)');
% ylabel('Frequency');
% title('Distance 0~100 nm')
% 
% figure;
% h200=histogram(diff_200,'BinWidth',0.2,'Normalization','Probability');
% xlim([-6 -0.5]);
% xlabel('log(diffusion coefficient) (\mum^2/s)');
% ylabel('Frequency');
% title('Distance 100~200 nm')
% 
% figure;
% h300=histogram(diff_300,'BinWidth',0.2,'Normalization','Probability');
% xlim([-6 -0.5]);
% xlabel('log(diffusion coefficient) (\mum^2/s)');
% ylabel('Frequency');
% title('Distance 200~300 nm')
% 
% figure;
% h400=histogram(diff_400,'BinWidth',0.2,'Normalization','Probability');
% xlim([-6 -0.5]);
% xlabel('log(diffusion coefficient) (\mum^2/s)');
% ylabel('Frequency');
% title('Distance 300~400 nm')
% 
% figure;
% h500=histogram(diff_500,'BinWidth',0.2,'Normalization','Probability');
% xlim([-6 -0.5]);
% xlabel('log(diffusion coefficient) (\mum^2/s)');
% ylabel('Frequency');
% title('Distance 400~500 nm')
% 
% figure;
% h600=histogram(diff_600,'BinWidth',0.2,'Normalization','Probability');
% xlim([-6 -0.5]);
% xlabel('log(diffusion coefficient) (\mum^2/s)');
% ylabel('Frequency');
% title('Distance 500~600 nm')
% 
% figure;
% h700=histogram(diff_700,'BinWidth',0.2,'Normalization','Probability');
% xlim([-6 -0.5]);
% xlabel('log(diffusion coefficient) (\mum^2/s)');
% ylabel('Frequency');
% title('Distance 600~700 nm')
% 
% figure;
% h800=histogram(diff_800,'BinWidth',0.2,'Normalization','Probability');
% xlim([-6 -0.5]);
% xlabel('log(diffusion coefficient) (\mum^2/s)');
% ylabel('Frequency');
% title('Distance 700~800 nm')
% 
% figure;
% h900=histogram(diff_900,'BinWidth',0.2,'Normalization','Probability');
% xlim([-6 -0.5]);
% xlabel('log(diffusion coefficient) (\mum^2/s)');
% ylabel('Frequency');
% title('Distance 800~900 nm')
% 
% figure;
% h1000=histogram(diff_1000,'BinWidth',0.2,'Normalization','Probability');
% xlim([-6 -0.5]);
% xlabel('log(diffusion coefficient) (\mum^2/s)');
% ylabel('Frequency');
% title('Distance 900~1000 nm')
% 
% figure;
% h1100=histogram(diff_1100,'BinWidth',0.2,'Normalization','Probability');
% xlim([-6 -0.5]);
% xlabel('log(diffusion coefficient) (\mum^2/s)');
% ylabel('Frequency');
% title('Distance 1000~1100 nm')
% 
% figure;
% h1200=histogram(diff_1200,'BinWidth',0.2,'Normalization','Probability');
% xlim([-6 -0.5]);
% xlabel('log(diffusion coefficient) (\mum^2/s)');
% ylabel('Frequency');
% title('Distance 1100~1200 nm')
% 
% figure;
% h1300=histogram(diff_1300,'BinWidth',0.2,'Normalization','Probability');
% xlim([-6 -0.5]);
% xlabel('log(diffusion coefficient) (\mum^2/s)');
% ylabel('Frequency');
% title('Distance 1200~1300 nm')
% 
% figure;
% h1400=histogram(diff_1400,'BinWidth',0.2,'Normalization','Probability');
% xlim([-6 -0.5]);
% xlabel('log(diffusion coefficient) (\mum^2/s)');
% ylabel('Frequency');
% title('Distance 1300~1400 nm')
% 
% figure;
% h1500=histogram(diff_1500,'BinWidth',0.2,'Normalization','Probability');
% xlim([-6 -0.5]);
% xlabel('log(diffusion coefficient) (\mum^2/s)');
% ylabel('Frequency');
% title('Distance 1400~1500 nm')
% 
% figure;
% h1600=histogram(diff_1600,'BinWidth',0.2,'Normalization','Probability');
% xlim([-6 -0.5]);
% xlabel('log(diffusion coefficient) (\mum^2/s)');
% ylabel('Frequency');
% title('Distance 1500~1600 nm')
% 
% figure;
% h1700=histogram(diff_1700,'BinWidth',0.2,'Normalization','Probability');
% xlim([-6 -0.5]);
% xlabel('log(diffusion coefficient) (\mum^2/s)');
% ylabel('Frequency');
% title('Distance 1600~1700 nm')
% 
% figure;
% h1800=histogram(diff_1800,'BinWidth',0.2,'Normalization','Probability');
% xlim([-6 -0.5]);
% xlabel('log(diffusion coefficient) (\mum^2/s)');
% ylabel('Frequency');
% title('Distance 1700~1800 nm')
% 
% figure;
% h1900=histogram(diff_1900,'BinWidth',0.2,'Normalization','Probability');
% xlim([-6 -0.5]);
% xlabel('log(diffusion coefficient) (\mum^2/s)');
% ylabel('Frequency');
% title('Distance 1800~1900 nm')
% 
% figure;
% h2000=histogram(diff_2000,'BinWidth',0.2,'Normalization','Probability');
% xlim([-6 -0.5]);
% xlabel('log(diffusion coefficient) (\mum^2/s)');
% ylabel('Frequency');
% title('Distance 1900~2000 nm')
% 
% 
