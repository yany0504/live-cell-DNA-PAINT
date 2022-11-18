%calculate tauD and tauB for each synapse in the
%image. V2: ROI became a sphere.
% V4: compare pre-LTP and post-LTP


%% Step 1. loading files
clc; 
clear;
close all;
size_synapse=500; % radius (distance from Homer)
exposure=0.1; %in sec
dimension=2; % 2 for 2D 3 for 3D

tauD_R1_1=33.284*2; % calibration of single 5xR1
tauD_R3_1=43.069*2; % calibration of single 7xR3

FilePath = 'E:\Data\DNA-PAINT(Live-cell)\Counting_qPAINT\cLTP_counting\cLTP_Mg-free\210604_qPAINT_cLTP';
OutputFile='qPAINT';
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

if ~exist('FileName_AMPAR_pre','var')|| isempty(FileName_AMPAR_pre)
    [userfilein, userdirin]=uigetfile({
         '*.csv','Data file (*.csv)';...
        '*.*','All Files (*.*)'},'Select PAINT file of pre-LTP AMPAR to process',...
        FilePath);
    FileName_AMPAR_pre=fullfile(userdirin,userfilein);
else
    if ~exist(FileName_AMPAR_pre,'file')
        fprintf('File not found: %s\n',FileName_AMPAR_pre);
        return;
    end
end
PAINT_pre=xlsread(FileName_AMPAR_pre);

if ~exist('FileName_AMPAR_post','var')|| isempty(FileName_AMPAR_post)
    [userfilein, userdirin]=uigetfile({
         '*.csv','Data file (*.csv)';...
        '*.*','All Files (*.*)'},'Select PAINT file of post-LTP AMPAR to process',...
        FilePath);
    FileName_AMPAR_post=fullfile(userdirin,userfilein);
else
    if ~exist(FileName_AMPAR_post,'file')
        fprintf('File not found: %s\n',FileName_AMPAR_post);
        return;
    end
end
PAINT_post=xlsread(FileName_AMPAR_post);


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
events_synapse_pre=cell(length(center_ROI),1);
if dimension == 2;
    for i=1:length(center_ROI);
        for j=1:length(PAINT_pre);
            if sqrt((center_ROI(i,1)-PAINT_pre(j,3))^2+(center_ROI(i,2)-PAINT_pre(j,4))^2) < size_synapse;
               events_synapse_pre{i}(end+1,:)=PAINT_pre(j,:);
            else
            end
        end
    end
    
elseif dimension == 3;
        for i=1:length(center_ROI);
        for j=1:length(PAINT_pre);
            if sqrt((center_ROI(i,1)-PAINT_pre(j,3))^2+(center_ROI(i,2)-PAINT_pre(j,4))^2+(center_ROI(i,3)-PAINT_pre(j,5))^2) < size_synapse;
               events_synapse_pre{i}(end+1,:)=PAINT_pre(j,:);
            else
            end
        end
    end
else
    error('dimension is wrong');
end

events_synapse_post=cell(length(center_ROI),1);
if dimension == 2;
    for i=1:length(center_ROI);
        for j=1:length(PAINT_post);
            if sqrt((center_ROI(i,1)-PAINT_post(j,3))^2+(center_ROI(i,2)-PAINT_post(j,4))^2) < size_synapse;
               events_synapse_post{i}(end+1,:)=PAINT_post(j,:);
            else
            end
        end
    end
    
elseif dimension == 3;
        for i=1:length(center_ROI);
        for j=1:length(PAINT_post);
            if sqrt((center_ROI(i,1)-PAINT_post(j,3))^2+(center_ROI(i,2)-PAINT_post(j,4))^2+(center_ROI(i,3)-PAINT_post(j,5))^2) < size_synapse;
               events_synapse_post{i}(end+1,:)=PAINT_post(j,:);
            else
            end
        end
    end
else
    error('dimension is wrong');
end


%% Step 4. calculate tauD and tauB for each synapse.
DarkTime_pre=zeros(length(events_synapse_pre),3);
DarkTime_pre(:,3)=1:length(DarkTime_pre);
BrightTime_pre=zeros(length(events_synapse_pre),3);
BrightTime_pre(:,3)=1:length(BrightTime_pre);
for i=1:length(events_synapse_pre);
    if size(events_synapse_pre{i},1)<3 %isempty(events_synapse{i});
        fprintf('Not enough events for the synapse.\n')
    else
        [DarkTime_pre(i,1:2) BrightTime_pre(i,1:2)]=CalculateTauDark(events_synapse_pre{i},exposure);
    end
    
end
clear i j;
% DarkTime_pre=DarkTime_pre(DarkTime_pre(:,1)~=0,:);
% BrightTime_pre=BrightTime_pre(BrightTime_pre(:,1)~=0,:);

DarkTime_post=zeros(length(events_synapse_post),3);
DarkTime_post(:,3)=1:length(DarkTime_post);
BrightTime_post=zeros(length(events_synapse_post),3);
BrightTime_post(:,3)=1:length(BrightTime_post);
for i=1:length(events_synapse_post);
    if size(events_synapse_post{i},1)<3 %isempty(events_synapse{i});
        fprintf('Not enough events for the synapse.\n')
    else
        [DarkTime_post(i,1:2) BrightTime_post(i,1:2)]=CalculateTauDark(events_synapse_post{i},exposure);
    end
    
end
clear i j;
% DarkTime_post=DarkTime_post(DarkTime_post(:,1)~=0,:);
% BrightTime_post=BrightTime_post(BrightTime_post(:,1)~=0,:);

%% convert to number of molecules

Number_pre=[tauD_R1_1./DarkTime_pre(:,1),DarkTime_pre(:,3)];
Number_pre(Number_pre(:,1)<0.5,1)=0;
Number_post=[tauD_R3_1./DarkTime_post(:,1),DarkTime_post(:,3)];
Number_post(Number_post(:,1)<0.5,1)=0;
for i=1:length(Number_pre);
    if Number_pre(i,1)==Inf;
        Number_pre(i,1)=0;
    else
    end
end
for i=1:length(Number_post);
    if Number_post(i,1)==Inf;
        Number_post(i,1)=0;
    else
    end
end
% filter out smaller than 1?
% Number_pre_filtered=Number_pre(Number_pre(:,1)>=1,:);
% Number_post_filtered=Number_post(Number_post(:,1)>=1,:);

% need to make synapse ID; # pre; # post; ratio;  
Number_synapse=[Number_pre(:,2),Number_pre(:,1),Number_post(:,1),Number_post(:,1)./Number_pre(:,1)];
Number_total=Number_pre(:,1)+Number_post(:,1);

non_zero_ratio=Number_synapse(Number_synapse(:,4)~=0&Number_synapse(:,4)~=Inf&~isnan(Number_synapse(:,4)),4);
noresponse=sum(Number_synapse(:,4)==0);
silentsynapses=sum(Number_synapse(:,4)==Inf);
noresponse_silent=sum(isnan(Number_synapse(:,4)));
non_zero_ratio=non_zero_ratio(non_zero_ratio<100);
mean_ratio=mean(non_zero_ratio)
std_ratio=std(non_zero_ratio)

figure;
h=histogram(non_zero_ratio);
ax=gca;
ax.FontSize=16;
xlabel('Fold increase # GluA1','FontSize',20)
ylabel('Count','FontSize',20)
title(['Mean = ' num2str(mean_ratio) ' \pm ' num2str(std_ratio) ' fold increased']);
Number_total=Number_total(Number_total<50);
mean_total=mean(Number_total)
std_total=std(Number_total)
figure;
h2=histogram(Number_total);
ax=gca;
ax.FontSize=16;
xlabel('Total # GluA1','FontSize',20)
ylabel('Count','FontSize',20)
title(['Mean = ' num2str(mean_total) ' \pm ' num2str(std_total)]);

