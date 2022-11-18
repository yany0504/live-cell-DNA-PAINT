%%Image contour to get cell area from cell image e.g. FP filling, Homer,
%%Actin, etc...

%% Filepath
clc; clear; close all;
FilePath='E:\Data\DNA-PAINT(Live-cell)\SPT_tracking\Basal_tracking\210204_livecellDNA-PAINT_HP_antiHomer1-mGeos_7xR3-scFv_R3_9nt-LD655\6\';
ImageFile='6_Homer-mGeos.tif';
ThunderSTORM='6_scFv-GluA1_corrected.csv';
outputname='filtered';
data=xlsread([FilePath ThunderSTORM]);
pixelsize=106;
%% make contour for the cell
I=imread([FilePath ImageFile]);
BW=imbinarize(I,'adaptive','Sensitivity',0.7);
%BW=im2bw(I,0.012);
figure;
imshowpair(I,BW,'montage');
figure;
[C, h]= imcontour(BW,1); % in h.ZData, the area inside the contour are 1 outside 0.

%% filter ThunderSTORM data inside the contours
x1=data(:,3)./pixelsize;
y1=data(:,4)./pixelsize;
x_pixel=round(x1);
y_pixel=round(y1);
% figure;
% scatter(x_pixel,y_pixel,'.');
datasize=size(data);
data_filtered=zeros(1,datasize(1,2));
for i=1:length(x_pixel);
    if h.ZData(y_pixel(i),x_pixel(i))==1;
        data_filtered=[data_filtered; data(i,:)];
    else
    end
end
csvwrite(strcat([FilePath ThunderSTORM(1:end-4) '_' outputname '.csv']),data_filtered);
        