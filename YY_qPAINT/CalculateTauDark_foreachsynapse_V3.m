%calculate tauD and tauB for each synapse in the
%image. V2: ROI became a sphere.

%% Step 1. loading files
clc; 
clear;
close all;
%clearvars -except cell1 cell2;%close all;
%parameters
size_synapse=500; % radius (distance from Homer)
exposure=0.1; %in sec
dimension=2; % 2 for 2D 3 for 3D
filterframe=0; % filter frame 1 on 0 off.
framemin=1;
framemax=20000;

tauD_R1_1=1; % calibration of single 5xR1
tauD_R3_1=67.24; % calibration of single 7xR3

FilePath = 'E:\Data\DNA-PAINT(Live-cell)\Counting_qPAINT\cLTP_counting\cLTP_Mg-free\210604_qPAINT_cLTP\1_cLTP\';
OutputFile='cLTP_tauD67.24s_qPAINT';
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

%% Step 1.5. filtering by frames
if filterframe==1;
    PAINT=PAINT(PAINT(:,2)>=framemin & PAINT(:,2)<=framemax,:);

else
end
%% Step 2. Setting ROI for each synapse
% center_ROI=zeros(length(Synapses1),3);

for i=1:length(Synapses1);
    center_x(i,1)=mean(Synapses1{i}(:,1));
    center_y(i,1)=mean(Synapses1{i}(:,2));
    center_z(i,1)=mean(Synapses1{i}(:,3));
end
center_ROI=[center_x center_y center_z];
clear i;
figure
for n=1:length(Synapses1)
        scatter3(Synapses1{n}(:,1),Synapses1{n}(:,2),Synapses1{n}(:,3),20,'filled');
        a=sprintf('%d',n);
        text(center_ROI(n,1)+100,center_ROI(n,2),center_ROI(n,3),a);
        hold all
end
hold off

axis equal
view(-360,-90);
%% Step 3. grouping PAINT results for each ROI.
events_synapse=cell(length(center_ROI),1);
if dimension == 2;
    for i=1:length(center_ROI);
        for j=1:length(PAINT);
            if sqrt((center_ROI(i,1)-PAINT(j,3))^2+(center_ROI(i,2)-PAINT(j,4))^2) < size_synapse;
               events_synapse{i}(end+1,:)=PAINT(j,:);
            else
            end
        end
    end
    
elseif dimension == 3;
        for i=1:length(center_ROI);
        for j=1:length(PAINT);
            if sqrt((center_ROI(i,1)-PAINT(j,3))^2+(center_ROI(i,2)-PAINT(j,4))^2+(center_ROI(i,3)-PAINT(j,5))^2) < size_synapse;
               events_synapse{i}(end+1,:)=PAINT(j,:);
            else
            end
        end
    end
else
    error('dimension is wrong');
end
%% cluster analysis for each synapse
list_synapse=[36, 37, 43, 63,45,42,49,7] ; % synapse number
%list_synapse=[43];
epsilon = 10;
minpts = 25;
darktime=cell(length(list_synapse),1);
for k=1:length(list_synapse);
    i_synapse=list_synapse(k);
    data=events_synapse{i_synapse};
    idx=dbscan(data(:,3:4),epsilon,minpts);
    data_temp=data;
    for i=1:length(data);
        if idx(i)==-1;
            data_temp(i,:)=zeros(1,9);
        else
        end
    end
    data_filtered=data_temp(data_temp(:,9)~=0,:);
    idx_filtered=idx(idx~=-1);
    figure;
    ax(1)=subplot(1,2,1);
    gscatter(data_filtered(:,3),data_filtered(:,4),idx_filtered);axis equal; hold on;
    legend('off')
    xlim([center_ROI(i_synapse,1)-500 center_ROI(i_synapse,1)+500])
    ylim([center_ROI(i_synapse,2)-500 center_ROI(i_synapse,2)+500])
    view(-360,-90);
    % calculate time for the synapse
    event_time=cell(max(idx_filtered),1);
    for j=1:max(idx_filtered);
        event_time{j,1}=data_filtered(idx_filtered==j,:);
    end
    DarkTime=zeros(length(event_time),3);
    DarkTime(:,3)=[1:length(DarkTime)]';
    BrightTime=zeros(length(event_time),3);
    BrightTime(:,3)=1:length(BrightTime);
    for i=1:length(event_time);
        if size(event_time{i},1)<5 %isempty(events_synapse{i});
            fprintf('Not enough events for the synapse.\n')
        else
            [DarkTime(i,1:2) BrightTime(i,1:2)]=CalculateTauDark(event_time{i},exposure);
        end
  
    end
    DarkTime=DarkTime(DarkTime(:,1)~=0,:);
    BrightTime=BrightTime(BrightTime(:,1)~=0,:);
    darktime{k}=DarkTime;
    ax(2)=subplot(1,2,2);
    scatplot(data_filtered(:,3),data_filtered(:,4));axis equal; 
    xlim([center_ROI(i_synapse,1)-500 center_ROI(i_synapse,1)+500])
    ylim([center_ROI(i_synapse,2)-500 center_ROI(i_synapse,2)+500])
    view(-360,-90);
    linkaxes(ax, 'xy');

end

darktime_sum=[];

for ii=1:length(list_synapse);
    darktime_sum=[darktime_sum; darktime{ii}];
end


mean(darktime_sum(:,1))
std(darktime_sum(:,1))

%%
figure;
scatplot(data(:,3),data(:,4));




