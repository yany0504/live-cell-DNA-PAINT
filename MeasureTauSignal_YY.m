% measuring tau fluorescence intensity in the proximity of Homer
% localizations. For better result, using background subtraction in ImageJ
% for tau image is helpful.
%% load files
clear; close all; clc;
%range=500; % diameter of ROI around Homer localization
pixelsize=107; % in nm.
mask=[0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,1,1,1,1,1,0,0,0;
    0,0,1,1,1,1,1,1,1,0,0;
    0,1,1,1,1,1,1,1,1,1,0;
    0,1,1,1,1,1,1,1,1,1,0;
    0,1,1,1,1,1,1,1,1,1,0;
    0,1,1,1,1,1,1,1,1,1,0;
    0,1,1,1,1,1,1,1,1,1,0;
    0,0,1,1,1,1,1,1,1,0,0;
    0,0,0,1,1,1,1,1,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0;];
path = "E:\Data\DNA-PAINT(Live-cell)\SPT_tracking\cLTP_tracking\210514_cLTP_HP_antiHomer1-EYFP\3_mGeos_cLTP\Cell1";
if ~exist('Synapse','var')|| isempty(Synapse)
    [userfilein1, userdirin1]=uigetfile({
         '*.mat','Data file (*.mat)';...
        '*.*','All Files (*.*)'},'Select Synapses1 file to process',...
        path);
    Synapse=fullfile(userdirin1,userfilein1);
else
    if ~exist(Synapse,'file')
        fprintf('File not found: %s\n',Synapse);
        return;
    else [userdirin1,~,~]=fileparts(Synapse);
        userdirin1=strcat(userdirin1,'\');
    end
end

if ~exist('tau','var')|| isempty(tau)
    [userfilein2, userdirin2]=uigetfile({
         '*.tif','Image file (*.tif)';...
        '*.*','All Files (*.*)'},'Select tau image to process',...
        path);
    tau=fullfile(userdirin2,userfilein2);
else
    if ~exist(tau,'file')
        fprintf('File not found: %s\n',tau);
        return;
    else [userdirin2,~,~]=fileparts(tau);
        userdirin2=strcat(userdirin2,'\');
    end
end
load([userdirin1 userfilein1]);
tau_IMG=imread([userdirin2 userfilein2]);
tau_IMG=double(tau_IMG);
figure(1);imshow(tau_IMG,[]);
figure(2);imshow(mask);
%% background subtraction

%% 2D center of synapses
center=zeros(length(Synapses1),2);
for i=1:length(Synapses1);
    center(i,:)=[mean(Synapses1{i,1}(:,1)),mean(Synapses1{i,1}(:,2))];
end
center_px=center/pixelsize;
center_px_ceil=ceil(center_px);
figure;imshow(tau_IMG,[]); hold on; scatter(center_px_ceil(:,1),center_px_ceil(:,2));
%% measure the intensity around the centers
tau_intensity=zeros(length(Synapses1),1);
for i=1:length(center_px_ceil);
    ROI=tau_IMG(center_px_ceil(i,1)-5:center_px_ceil(i,1)+5,center_px_ceil(i,2)-5:center_px_ceil(i,2)+5);
    intensity=ROI.*mask;
    tau_intensity(i,1)=sum(intensity,'all')/sum(mask,'all');
end

