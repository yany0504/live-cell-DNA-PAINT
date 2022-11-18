% All steps analysis for DNA-PAINT. You don't do
% tracking but do cluster analysis and distance
% measurement.
%% Step 0 Parameters. YOU NEED TO TAKE A LOOK AT THIS PART CAREFULLY, OTHERWISE YOU WILL GET ERRORS OR MISINTERPRETED RESULTS.
clc;clear;close all;% Close all figures

path = 'H:\Data_from07222021\DNA-PAINT(Live-cell)\Counting_qPAINT\CTRL_counting\201106_Live-cell_DNA-PAINT_precisionofcounting_fixedHP\1\';
outputpath='H:\Data_from07222021\DNA-PAINT(Live-cell)\Counting_qPAINT\CTRL_counting\201106_Live-cell_DNA-PAINT_precisionofcounting_fixedHP\1\';
pathsplit=strsplit(path,'\');
outputindex=pathsplit(end-1);
outputname='1_distanceonly'; % this will go to all the output filex e.g. result_tracking_picasso.txt
isitpicasso=0;
actualdimension=2;
% frame selection
timerange=0; %if you want to specify frame = 1, otherwise = 0.
mint=1;
maxt=3000;
% ROI setting in nm.
useROI=0; % if you want to specify ROI 1, otherwise 0.
xmin=5606.8;
xmax=13653.2;
ymin=0;
ymax=18404;
% parameters for z threshold and factor
threshold=600; %YY corrected for ThunderSTORM. originally 600 from line ~59
if isitpicasso==1;
    factor=1;
else
    factor=0.79; % for picasso if you have already multiplied the factor in localization step, this should be 1 for ThunderSTORM =0.79
end
% params for cluster analysis
k_cluster=50;
Eps_cluster=300;

%% Step1. load files

if ~exist('fileName1','var')|| isempty(fileName1)
    [userfilein, userdirin]=uigetfile({
         '*.csv','Data file (*.csv)';...
        '*.*','All Files (*.*)'},'Select the Homer file to process',...
        path);
    fileName1=fullfile(userdirin,userfilein);
else
    if ~exist(fileName1,'file')
        fprintf('File not found: %s\n',fileName1);
        return;
    end
end
if ~exist('fileName2','var')|| isempty(fileName2)
    [userfilein, userdirin]=uigetfile({
         '*.csv','Data file (*.csv)';...
        '*.*','All Files (*.*)'},'Select the DNA-PAINT file to process',...
        path);
    fileName2=fullfile(userdirin,userfilein);
else
    if ~exist(fileName2,'file')
        fprintf('File not found: %s\n',fileName2);
        return;
    else [userdirin,~,~]=fileparts(fileName2);
        userdirin=strcat(userdirin,'\');
    end
end

if ~exist('fileName3','var')|| isempty(fileName3)
    [userfilein, userdirin]=uigetfile({
         '*.csv','Data file (*.csv)';...
        '*.*','All Files (*.*)'},'Select the Bassoon file to process',...
        path);
    fileName3=fullfile(userdirin,userfilein);
else
    if ~exist(fileName3,'file')
        fprintf('File not found: %s\n',fileName3);
        return;
    end
end
presynapse=xlsread(fileName3);
% if ~exist('threshold','var')||isempty(threshold)||~exist('factor','var')||isempty(factor)
%     prompt = {'Enter threshold:','Enter factor of z:'};
%     dlg_title = 'Input for preparing data';
%     num_lines = 1;
%     defaultanswer={'600','0.79'};
%     answer = inputdlg(prompt,dlg_title,num_lines,defaultanswer);
%     threshold=str2double(answer(1));
%     factor=str2double(answer(2));
% end
% threshold=1000; %YY corrected for ThunderSTORM. originally 600
% factor=0.79;

data_PALM=xlsread(fileName1);
data_PALM(isnan(data_PALM))=Inf;
%%% Anne's revision
if actualdimension==3;
    data_PALM=data_PALM(abs(data_PALM(:,5))<threshold,:); %for 3D
    data_PALM(:,5)=data_PALM(:,5)*factor; %for 3D
else actualdimension==2;
end

%%%%%%%%%%%%%% revision
fid=fopen([userdirin 'data_PALM_' outputname '.txt'],'w');
fprintf(fid,'%d %f %f %f %f %f %f\n',[data_PALM(:,1) data_PALM(:,2:7)]');    
fclose(fid);
%%%%%%%%%%%%%% revision
data_QD=xlsread(fileName2);
data_QD(isnan(data_QD))=Inf;
if actualdimension==3;
    data_QD=data_QD(abs(data_QD(:,5))<threshold,:);
    data_QD(:,5)=data_QD(:,5)*factor;
else actualdimension==2;
end
data_QD=sortrows(data_QD,2);


if timerange ==1;
    idx_t=(data_QD(:,2)>=mint)&(data_QD(:,2)<=maxt);
    indices=find(idx);
    for i=(1:length(indices))
         data_QD_temp(i,:)=data_QD(indices(i),:);
    end
    data_QD=data_QD_temp;
else
end



% data filtering for ROI.
if useROI == 1;
    idx_QD=find(data_QD(:,3) > xmin & data_QD(:,3) < xmax & data_QD(:,4) > ymin & data_QD(:,4) < ymax);
    idx_PALM=find(data_PALM(:,3) > xmin & data_PALM(:,3) < xmax & data_PALM(:,4) > ymin & data_PALM(:,4) < ymax);
    QD=[];
    for i=(1:length(idx_QD))
        QD(i,:)=data_QD(idx_QD(i),:);
    end
    PALM=[];
    for i=(1:length(idx_PALM))
        PALM(i,:)=data_PALM(idx_PALM(i),:);
    end
    data_QD=QD;
    data_PALM=PALM;
else
end

%% Step2. cluster analysis for defining synapses.

if actualdimension==3;
    palmX=data_PALM(:,3);
    palmY=data_PALM(:,4);
    palmZ=data_PALM(:,5);
    data_Syn=[palmX, palmY, palmZ];
else actualdimension==2;
    palmX=data_PALM(:,3);
    palmY=data_PALM(:,4);
    palmZ=zeros(length(palmX),1);
    data_Syn=[palmX, palmY, palmZ];
end
    
% Finding clusters
[Class,type]=dbscan_conservative(data_Syn,k_cluster,Eps_cluster); 
% Make new matrix  
Syn=[palmX,palmY,palmZ,Class',type'];
% Seperate and plot clusters
x=Syn(:,1);
y=Syn(:,2);
z=Syn(:,3); 
particle=Syn(:,4);
PtType=Syn(:,5);
d=length(particle); % total row number
ptotal=max(particle);% total cluster number
numArrays=ptotal;
Synapses1=cell(numArrays,1);
figure
for n=1:ptotal
        dp=find(particle==n); % matrix indcis of points in cluster #k 
        dplength=length(dp);
        %dpmax=max(dp);
        %dpmin=min(dp);
        m(n)=n; %test k value
        x1=x(dp);
        y1=y(dp);
        z1=z(dp);
        PtType1=PtType(dp);
        PtNum=n*ones(dplength,1);
        Syn1=[x1,y1,z1,PtType1,PtNum]; 
        Synapses1{n}=Syn1;   % Synapses: X, Y, Z, Type, Partical ID
        %eval(['Synapse_' num2str(k) '=[x1,y1,z1,PtType1]']);% creating submatrix for particle# k sub_k 
        scatter3(x1,y1,z1,20,'filled');
        hold all
        clear dp dplength dpmax dpmin x1 y1 z1 PtType1
end

hold off
clear x y z PtType particle k d m x1 y1 z1 dp Syn PtNum 
saveas(gcf,strcat(userdirin,['syn_' outputname]),'fig');
save(strcat(userdirin,['Synapses1_' outputname '.mat']),'Synapses1');
axis equal


%% Step3. Measure distance between PAINT and Synapses.
if actualdimension==3;
    x=data_QD(:,3);
    y=data_QD(:,4);
    z=data_QD(:,5);
elseif actualdimension==2;
    x=data_QD(:,3);
    y=data_QD(:,4);
    z=zeros(length(data_QD),1);
else
    printf('Dimension is wrong');
end

    for n=1:length(x)
    position(n,1)=x(n);
    position(n,2)=y(n);
    position(n,3)=z(n);
    nSyn=length(Synapses1);
    for rg=1:nSyn
        syn_rg=Synapses1{rg}; % load data
        % Find center of the synapse
        clear syn_position
        Ctr0=mean(syn_rg);; % Find center of a subcluster and the range of influence of the center
        %dist=zeros(length(traces{1}));
        % Find traces close to the synapse
        syn_position(rg,1)=(Ctr0(1)-position(n,1)).*(Ctr0(1)-position(n,1));
        syn_position(rg,2)=(Ctr0(2)-position(n,2)).*(Ctr0(2)-position(n,2));
        syn_position(rg,3)=(Ctr0(3)-position(n,3)).*(Ctr0(3)-position(n,3));
        distance(rg)=sqrt(syn_position(rg,1)+syn_position(rg,2)+syn_position(rg,3))/1000;
    end
dist = [distance];
dist2 = dist';
mini=min(dist2);
t= find(dist2 == mini);
a=Synapses1{t};
syn_cent0= mean(a);
Syn_ctrX = syn_cent0(1);
Syn_ctrY = syn_cent0(2);
Syn_ctrZ = syn_cent0(3);
TS_X=(position(n,1)-Syn_ctrX)/1000; % Trace to synapse distance delta X
TS_Y=(position(n,2)-Syn_ctrY)/1000; 
TS_Z=(position(n,3)-Syn_ctrZ)/1000;
TS_dist{n}=sqrt(TS_X.*TS_X+TS_Y.*TS_Y+TS_Z.*TS_Z);
end

Dist_measure=cell2mat(TS_dist');
DistBin=(0:0.02:2);
figure;
subplot(1,2,1);
hDistance=histogram(Dist_measure,DistBin,'Normalization','pdf');
xlim([0 2]);
xlabel('Distance (\mum), PDF','FontSize',18);
ylabel('Frequency','FontSize',18);
subplot(1,2,2);
hDistanceCDF=histogram(Dist_measure,DistBin,'Normalization','cdf');
xlim([0 2]);
xlabel('Distance (\mum), CDF','FontSize',18);
ylabel('CDF','FontSize',18);
set(gcf, 'Position', [100 100 1600 500]);
saveas(gcf,strcat(outputpath,['TS_dist_' outputname]),'fig');
xlswrite(strcat(outputpath,['TS_dist_cal_' [outputname '_' outputindex{1}] '.xlsx']), Dist_measure);

%% cluster analysis for Bassoon.

data_pre=[presynapse(:,3:4) zeros(length(presynapse),1)];
% Finding clusters
[Class,type]=dbscan_conservative(data_pre,k_cluster,Eps_cluster); 
% Make new matrix  
Syn=[data_pre(:,1),data_pre(:,2),data_pre(:,3),Class',type'];
% Seperate and plot clusters
x=Syn(:,1);
y=Syn(:,2);
z=Syn(:,3); 
particle=Syn(:,4);
PtType=Syn(:,5);
d=length(particle); % total row number
ptotal=max(particle);% total cluster number
numArrays=ptotal;
PreSynapses1=cell(numArrays,1);
figure
for n=1:ptotal
        dp=find(particle==n); % matrix indcis of points in cluster #k 
        dplength=length(dp);
        %dpmax=max(dp);
        %dpmin=min(dp);
        m(n)=n; %test k value
        x1=x(dp);
        y1=y(dp);
        z1=z(dp);
        PtType1=PtType(dp);
        PtNum=n*ones(dplength,1);
        PreSyn1=[x1,y1,z1,PtType1,PtNum]; 
        PreSynapses1{n}=PreSyn1;   % Synapses: X, Y, Z, Type, Partical ID
        %eval(['Synapse_' num2str(k) '=[x1,y1,z1,PtType1]']);% creating submatrix for particle# k sub_k 
        scatter3(x1,y1,z1,20,'filled');
        hold all
        clear dp dplength dpmax dpmin x1 y1 z1 PtType1
end

hold off
clear x y z PtType particle k d m x1 y1 z1 dp Syn PtNum 
saveas(gcf,strcat(userdirin,['syn_' outputname]),'fig');
save(strcat(userdirin,['PreSynapses1_' outputname '.mat']),'PreSynapses1');
axis equal
%% Cetegorize Homer clusters into 1) active: which have Bassoon cluster within xxx nm and 2) inactive: which don't 
% label 1 for active and 2 for inactive.
centerSyn=zeros(length(Synapses1),3);
for i=1:length(Synapses1);
    centerSyn(i,:)=[mean(Synapses1{i,1}(:,1)),mean(Synapses1{i,1}(:,2)),mean(Synapses1{i,1}(:,3))];
end

centerpre=zeros(length(PreSynapses1),3);
for i=1:length(PreSynapses1);
    centerpre(i,:)=[mean(PreSynapses1{i,1}(:,1)),mean(PreSynapses1{i,1}(:,2)),mean(PreSynapses1{i,1}(:,3))];
end

minDist=zeros(length(centerSyn),1);
for i=1:length(minDist);
    dist2=((centerpre-centerSyn(i,:)).^2);    
    dist=sqrt(dist2(:,1)+dist2(:,2)+dist2(:,3));
    minDist(i)=min(dist);
end

save([outputpath outputname '_Distance_Homer_Bassoon.mat'],'minDist');



