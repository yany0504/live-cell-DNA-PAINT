%% loading files
clear;clc;
FilePath ='E:\Data\DNA-Tracking\201016_Live-DNA-PAINT_Tracking\4_actualDatawith8nt\'
ThunderSTORM='AMPAR_All_c.csv'
data=xlsread([FilePath ThunderSTORM]);
outputtitle='8-nt imager'

%% Histogram photon number
figure;
h_ph=histogram(data(:,8),'BinWidth',20);
xlabel('Photon number','FontSize',22);
ylabel('Count','FontSize',22);
xlim([0 5000]);
title([outputtitle ', Photon number'],'FontSize',18);
ax=gca;
ax.YAxis.FontSize=14;
ax.XAxis.FontSize=14;

%% Histogram BG noise
figure;
h_BG=histogram(data(:,10),'BinWidth',1);
xlabel('Background Noise (photon)','FontSize',22);
ylabel('Count','FontSize',22);
xlim([0 50]);
title([outputtitle ', Background Noise'],'FontSize',18);
ax=gca;
ax.YAxis.FontSize=14;
ax.XAxis.FontSize=14;

%% Histogram accuracy (lateral)
figure;
h_latacc=histogram(data(:,11),'BinWidth',1);
xlabel('Lateral Accuracy (nm)','FontSize',22);
ylabel('Count','FontSize',22);
xlim([0 50]);
title([outputtitle ', Lateral Accuracy'],'FontSize',18);
ax=gca;
ax.YAxis.FontSize=14;
ax.XAxis.FontSize=14;


%% Histogram accuracy (lateral)
figure;
h_axiacc=histogram(data(:,12),'BinWidth',1);
xlabel('Axial Accuracy (nm)','FontSize',22);
ylabel('Count','FontSize',22);
xlim([0 150]);
title([outputtitle ', Axial Accuracy'],'FontSize',18);
ax=gca;
ax.YAxis.FontSize=14;
ax.XAxis.FontSize=14;
