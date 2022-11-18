%calculate tauD and tauB for each synapse in the
%image.

%% Step 1. loading files
clc;clear;%close all;
%parameters
size_synapse=1000; % radius (distance from Homer)
exposure=0.1; %in sec
dimension=2 % 2 for 2D 3 for 3D

FilePath = 'E:\Data\DNA-PAINT\200817_DNAPAINT_PCAPCD_HP\';
OutputFile='IMG';
if ~exist('FileName_Synapse','var')|| isempty(FileName_Synapse)
    [userfilein, userdirin]=uigetfile({
         '*.mat','Data file (*.mat)';...
        '*.*','All Files (*.*)'},'Select Synapses1 file to process',...
        FilePath);
    FileName_Synapse=fullfile(userdirin,userfilein);
else
    if ~exist(FileName_Synapse,'file')
        fprintf('File not found: %s\n',FileName_Synapse);
        return;
    end
end
load(FileName_Synapse);
if ~exist('FileName_AMPAR','var')|| isempty(FileName_PAINT)
    [userfilein, userdirin]=uigetfile({
         '*.csv','Data file (*.csv)';...
        '*.*','All Files (*.*)'},'Select PAINT file to process',...
        FilePath);
    FileName_PAINT=fullfile(userdirin,userfilein);
else
    if ~exist(FileName_PAINT,'file')
        fprintf('File not found: %s\n',FileName_PAINT);
        return;
    end
end

PAINT=xlsread(FileName_PAINT);


%% Step 2. Setting ROI for each synapse
% center_ROI=zeros(length(Synapses1),3);
for i=1:length(Synapses1);
    center_x(i,1)=mean(Synapses1{i}(:,1));
    center_y(i,1)=mean(Synapses1{i}(:,2));
    center_z(i,1)=mean(Synapses1{i}(:,3));
end
center_ROI=[center_x center_y center_z];
clear i;
%x_min, x_max, y_min, y_max, z_min, z_max
ROI=[center_x-size_synapse center_x+size_synapse center_y-size_synapse center_y+size_synapse center_z-size_synapse center_z+size_synapse];
%% Step 3. grouping PAINT results for each ROI.
events_synapse=cell(length(ROI),1);
if dimension == 2;
    for i=1:length(ROI);
        for j=1:length(PAINT);
            if PAINT(j,3)>ROI(i,1) & PAINT(j,3)<ROI(i,2) & PAINT(j,4)>ROI(i,3) & PAINT(j,4)<ROI(i,4)
                events_synapse{i}(end+1,:)=PAINT(j,:);
            else
            end
        end
    end
    
elseif dimension == 3;
        for i=1:length(ROI);
        for j=1:length(PAINT);
            if PAINT(j,3)>ROI(i,1) & PAINT(j,3)<ROI(i,2) & PAINT(j,4)>ROI(i,3) & PAINT(j,4)<ROI(i,4) & PAINT(j,5)>ROI(i,5) & PAINT(j,5)<ROI(i,6)
                events_synapse{i}(end+1,:)=PAINT(j,:);
            else
            end
        end
    end
else
    error('dimension is wrong');
end

%% Step 4. calculate tauD and tauB for each synapse.
DarkTime=zeros(length(events_synapse),2);
BrightTime=zeros(length(events_synapse),2);
for i=1:length(events_synapse);
    [DarkTime(i,:) BrightTime(i,:)]=CalculateTauDark(events_synapse{i},exposure);
end

