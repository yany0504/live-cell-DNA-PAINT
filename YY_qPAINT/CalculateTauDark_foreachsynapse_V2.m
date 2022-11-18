%calculate tauD and tauB for each synapse in the
%image. V2: ROI became a sphere.

%% Step 1. loading files
clc; 
clear;
close all;
size_synapse=500; % radius (distance from Homer)
exposure=0.1; %in sec
dimension=2; % 2 for 2D 3 for 3D
filterframe=0; % filter frame 1 on 0 off.
framemin=1;
framemax=20000;

tauDark=77.5; % darktime of single puncta

FilePath = 'E:\Data\DNA-PAINT(Live-cell)\Counting_qPAINT\cLTP_counting\cLTP_singlelabeling\210623_qPAINT_phalloidin_cLTP_vehicle\cLTP\4\';
OutputFile='4_cLTP_qPAINT';
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

if ~exist('FileName_minDist','var')|| isempty(FileName_minDist)
    [userfilein, userdirin]=uigetfile({
         '*.mat','Data file (*.mat)';...
        '*.*','All Files (*.*)'},'Select Homer-Bassoon Distance file to process',...
        FilePath);
    FileName_minDist=fullfile(userdirin,userfilein);
else
    if ~exist(FileName_minDist,'file')
        fprintf('File not found: %s\n',FileName_minDist);
        return;
    end
end
load(FileName_minDist);

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

%% Step 4. calculate tauD and tauB for each synapse.
DarkTime=zeros(length(events_synapse),2);
DarkTimeCI=zeros(length(events_synapse),2);
DarkTime(:,2)=1:length(DarkTime);
BrightTime=zeros(length(events_synapse),2);
BrightTimeCI=zeros(length(events_synapse),2);
BrightTime(:,2)=1:length(BrightTime);


for i=1:length(events_synapse);
    if size(events_synapse{i},1)<3 %isempty(events_synapse{i});
        fprintf('Not enough events for the synapse.\n')
    else
        [DarkTime(i,1), DarkTimeCI(i,:), BrightTime(i,1), BrightTimeCI(i,:)]=CalculateTauDark(events_synapse{i},exposure,5);
    end
    
end
clear i j;
DarkTime=DarkTime(DarkTime(:,1)~=0,:);
DarkTimeCI=DarkTimeCI(DarkTime(:,1)~=0,:);
BrightTime=BrightTime(BrightTime(:,1)~=0,:);
BrightTimeCI=BrightTimeCI(BrightTime(:,1)~=0,:);
mean(DarkTime(:,1))
std(DarkTime(:,1))

%% Bonus step. Estimate the number of target in ROI using qPAINT.

Number=tauDark./DarkTime(:,1);


Number_mean=mean(Number(:,1))
Number_std=std(Number(:,1))
Number_SEM=Number_std/sqrt(length(Number(:,1)))



%% categorizing type of synapses 1) active = Homer-Bassoon distance < 250 nm, 2) inactive = Homer-Bassoon distance > 250 nm.

synapse_cat=double(minDist<250);
Number_active=Number(synapse_cat==1);
Number_inactive=Number(synapse_cat==0);
Number_active_mean=mean(Number_active)
Number_active_std=std(Number_active)
Number_active_SEM=Number_active_std/sqrt(length(Number_active(:,1)))
Number_inactive_mean=mean(Number_inactive)
Number_inactive_std=std(Number_inactive)
Number_inactive_SEM=Number_inactive_std/sqrt(length(Number_inactive(:,1)))

save([userdirin 'qPAINT_' OutputFile '.mat']);





