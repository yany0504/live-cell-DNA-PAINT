% plot intensity trace of blinking image in ROI based on the drift file
clear;close all;clc;

Path='H:\Data_from07222021\DNA-PAINT(Live-cell)\kinetics_calibration\210726_DNA-PAINT_calib_2x5x16xR1_bufferC_R1_7ntCy3B_500pM\2xR1\';
Image='2xR1_001.tif';
Drift='Drift.csv';
interp_factor=4;
pixel_size_ir=44;
pixel_size_ccd=107;
box_size=4;
ImageFile=[Path Image];
DriftFile=[Path Drift];
IMG_info=imfinfo(ImageFile);
drift=xlsread(DriftFile);
Height=IMG_info.Height;
Width=IMG_info.Width;
Frame_num=length(IMG_info);
%% read image;
IMG(:,:,1)=imread(ImageFile,1,'Info',IMG_info);

for i=2:Frame_num;
   IMG(:,:,i)=imread(ImageFile,i,'Info',IMG_info);
   i
end
IMG=double(IMG);
%% correct angle between IR and EMCCD
New_angle_camera_direct(drift,pixel_size_ir,pixel_size_ccd,interp_factor);
drift_traj=drift_angle(:,1:2)/pixel_size_ccd;

%% visualize image and select a point from the first frame;
xi=[];
yi=[];
dipshow(IMG(:,:,1),hot,'lin');
diptruesize(400);
[xi,yi] = getpts;
roi_x=[xi-box_size/2,xi+box_size/c2];
roi_y=[yi-box_size/2,yi+box_size/2];
roi_xt=zeros(Frame_num,2,length(roi_x));
roi_yt=zeros(Frame_num,2,length(roi_y));
for i=1:length(roi_x);
    roi_xt_temp=ones(Frame_num,2).*roi_x(i,:);
    roi_yt_temp=ones(Frame_num,2).*roi_y(i,:);
    roi_xt(:,:,i)=roi_xt_temp(:,:)-drift_angle(:,1);
    roi_yt(:,:,i)=roi_yt_temp(:,:)-drift_angle(:,2);
end
roi_xt=ceil(roi_xt);
roi_yt=ceil(roi_yt);

%% 
intensity_roi=zeros(Frame_num,length(roi_x));
for i=1:length(roi_x);
    IMG_roi=IMG(roi_xt(:,1,i):roi_xt(:,2,i)-1,roi_yt(:,1,i):roi_yt(:,2,i)-1,:);    
    intensity(:,i)=mean(mean(IMG_roi,2),1);
end

%% plot
close all
j=34;
x=1:Frame_num;
figure;plot(x,intensity(:,j));
dist_inten=intensity(:,j);
dist_inten=dist_inten(dist_inten>200);
figure;h1=histogram(dist_inten);
%% autocorrelation
LagTime=1000;
figure;
[acf,lags,bounds,h]=autocorr(intensity(:,j),'NumLags',LagTime,'NumSTD',5);

xlim([0.5 LagTime])
hold on;

f1=fit(lags,acf,'exp1');

h1=plot(f1,'k');
set(gca,'Xscale','log');

