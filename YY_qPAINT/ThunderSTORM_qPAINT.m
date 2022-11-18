%ThunderSTORM image visualizer.

%% loading files
clear;clc;close all;
FilePath ='H:\Data_from07222021\DNA-PAINT(Live-cell)\SPT_tracking\Basal_tracking\220217_live-PAINT_basal_37C_9nt\Cell3\';
ThunderSTORM='Cell3_GluA1_cycle1_corrected.csv';
%ThunderSTORM='GluA1_000.csv';
data=xlsread([FilePath ThunderSTORM]);
data_continue=[];
data_continue=csvread([FilePath ThunderSTORM],length(data)+1,0);
data_whole=[data;data_continue];
data_whole=data;
outputtitle='test';
height=512;
width=512;
CCD_pixelSize=107;
Magnification=5;
exposure=0.1; % second
timestep=1000; %frames
max_frame=25000;
data1=data_whole;
data1(:,3:4)=data1(:,3:4)./CCD_pixelSize;
data1=data1(data1(:,3)>0&data1(:,4)>0&data1(:,3)<width&data1(:,4)<height,:);

locIm = Loc2Img(data1(:,3:4),Magnification,height,width);
dipshow(locIm,hot,'lin');

roi_list=cell(100,1);
again = true;
regionCount = 0;

%%
while again && regionCount < 100
    promptMessage = sprintf('Draw region #%d in the upper right image,\nor Quit?', regionCount + 1);
	titleBarCaption = 'Continue?';
	button = questdlg(promptMessage, titleBarCaption, 'Draw', 'Quit', 'Draw');
	if strcmpi(button, 'Quit')
		break;
    end
    regionCount = regionCount + 1;
    message = sprintf(['Find ROI using pan (Ctrl+p) and zoom (Ctrl+z).\n If you are ready, click Ok, then draw ROI # ' num2str(regionCount)]);
	uiwait(msgbox(message));
    roi=images.roi.Polygon;
    draw(roi)
    roi_list{regionCount}=roi;
    
end
roi_list=roi_list(~cellfun('isempty',roi_list));

%% do qPAINT for each ROI.
DarkTime=zeros(length(roi_list),1);
DarkTimeCI=zeros(length(roi_list),2);
BrightTime=zeros(length(roi_list),1);
BrightTimeCI=zeros(length(roi_list),2);

for i=1:length(roi_list);
    in_roi=inpolygon(data1(:,3),data1(:,4),roi_list{i}.Position(:,1)/Magnification,roi_list{i}.Position(:,2)/Magnification);
    %sum(in_roi)
    data_ROI=data1(in_roi,:);
    [DarkTime(i), DarkTimeCI(i,:), BrightTime(i), BrightTimeCI(i,:)]=CalculateTauDark(data_ROI,exposure,5);
end

%% do qPAINT for each time window;
num_step=max_frame/timestep;
DarkTime_t=zeros(length(roi_list),1,num_step);
DarkTimeCI_t=zeros(length(roi_list),2,num_step);
BrightTime_t=zeros(length(roi_list),1,num_step);
BrightTimeCI_t=zeros(length(roi_list),2,num_step);

for i=1:num_step;
    data_time=data1(data1(:,2)>timestep*(i-1)&data1(:,2)<=timestep*i,:);
    for j=1:length(roi_list);
        in_roi=inpolygon(data_time(:,3),data_time(:,4),roi_list{j}.Position(:,1)/Magnification,roi_list{j}.Position(:,2)/Magnification);
        if sum(in_roi) == 0;
            DarkTime_t(j,1,i)=0;
            DarkTimeCI_t(j,:,i)=[0,0];
            BrightTime_t(j,1,i)=0;
            BrightTimeCI_t(j,:,i)=[0,0];
        else
            data_time_ROI=data_time(in_roi,:);
            [DarkTime_t(j,1,i), DarkTimeCI_t(j,:,i), BrightTime_t(j,1,i), BrightTimeCI_t(j,:,i)]=CalculateTauDark(data_time_ROI,exposure,5);
        end
    end
end
%% plot
DarkTime1=69.4;
number=DarkTime1./DarkTime_t;
numberCI=DarkTime1./DarkTimeCI_t;
number_SEM=numberCI(:,1,:)-number(:,1,:);

number_norm=number(:,:,:)./number(:,:,1);
number_SEM_norm=number_SEM./number(:,:,1);
i=3;
x=1:num_step;
figure;
errorbar(x,reshape(number_norm(i,1,:),[num_step,1]),reshape(number_SEM_norm(i,1,:),[num_step,1]));
hold on
scatter(x,reshape(number_norm(i,1,:),[num_step,1]));

% bulk assay;

[N, edges,bin]=histcounts(data_whole(:,2),'BinWidth',timestep);
