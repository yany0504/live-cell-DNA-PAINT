%% loading files
clear;clc;close all;
FilePath ='H:\Data_from07222021\DNA-PAINT(Live-cell)\kinetics_calibration\210426_DNA-PAINT_calib_ACSF_1x5x_16xR1_R1_8nt\16x1\';
ThunderSTORM='16xR1_1nM_R1_8nt-Cy3B_003_locs.csv';
%ThunderSTORM='GluA1_000.csv';
data=xlsread([FilePath ThunderSTORM]);
data_continue=[];
data_continue=csvread([FilePath ThunderSTORM],length(data)+1,0);
data_whole=[data;data_continue];
outputtitle='test';
height=256;
width=256;
CCD_pixelSize=107;
Magnification=5;
box_size=200;
exposure=0.1; % second
% timestep=1000; %frames
% max_frame=25000;
data1=data_whole;
data1(:,3:4)=data1(:,3:4)./CCD_pixelSize;
data1=data1(data1(:,3)>0&data1(:,4)>0&data1(:,3)<width&data1(:,4)<height,:);

locIm = Loc2Img(data1(:,3:4),Magnification,height,width);
dipshow(locIm,hot,'lin');

roi_list=cell(100,1);
again = true;
regionCount = 0;

%% visualize image and select a point from the first frame;
xi=[];
yi=[];
dipshow(locIm,hot,'lin');
diptruesize(400);

%% automatic peak finding

center=FastPeakFind(locIm,50);
xi=zeros(length(center)/2,1);
yi=zeros(length(center)/2,1);
for i=1:length(center)/2;
    xi(i,1)=center(2*i-1);
    yi(i,1)=center(2*i);
    
end
xi_nm=xi*CCD_pixelSize/Magnification;
yi_nm=yi*CCD_pixelSize/Magnification;
roi_x=[xi_nm-box_size/2,xi_nm+box_size/2];
roi_y=[yi_nm-box_size/2,yi_nm+box_size/2];
events_roi=cell(length(xi),1);

for i=1:length(xi);
    events_roi{i}=data_whole(data_whole(:,3)>roi_x(i,1)&data_whole(:,3)<roi_x(i,2)&data_whole(:,4)>roi_y(i,1)&data_whole(:,4)<roi_x(i,2),:);
end

events_roi=events_roi(~cellfun('isempty',events_roi));
%% manual peak finding
[xi,yi] = getpts;
xi_nm=xi*CCD_pixelSize/Magnification;
yi_nm=yi*CCD_pixelSize/Magnification;
roi_x=[xi_nm-box_size/2,xi_nm+box_size/2];
roi_y=[yi_nm-box_size/2,yi_nm+box_size/2];

events_roi=cell(length(xi),1);

for i=1:length(xi);
    events_roi{i}=data_whole(data_whole(:,3)>roi_x(i,1)&data_whole(:,3)<roi_x(i,2)&data_whole(:,4)>roi_y(i,1)&data_whole(:,4)<roi_x(i,2),:);
end

events_roi=events_roi(~cellfun('isempty',events_roi));

%% do qPAINT for each ROI.
DarkTime=zeros(length(events_roi),1);
DarkTimeCI=zeros(length(events_roi),2);
BrightTime=zeros(length(events_roi),1);
BrightTimeCI=zeros(length(events_roi),2);

for i=1:length(events_roi);

    [DarkTime(i), DarkTimeCI(i,:), BrightTime(i), BrightTimeCI(i,:)]=CalculateTauDark(events_roi{i},exposure,5);
end

figure;histogram(BrightTime);
figure;histogram(DarkTime);