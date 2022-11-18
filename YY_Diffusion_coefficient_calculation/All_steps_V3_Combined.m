function All_steps(fileName1,fileName2,threshold,factor,k,Eps,optional)
%line 39 threshold fixed 600 nm 0.79
% line 181 enter k and eps: fixed 50, 500
% Modified by Anne. Comapred to V2, V3 used tracking fucntion on mTraceIN
% and mTraceOUT in order to avoid connecting discontinued tracdes of same particle in or
% out of synapse
% Modified by Anne, 1/14/2014, added Dif_track ID, for note which trace ID is for corrisponding Dif_track

%% Step 0 Parameters. YOU NEED TO TAKE A LOOK AT THIS PART CAREFULLY, OTHERWISE YOU WILL GET ERRORS OR MISINTERPRETED RESULTS.
clc;clear;
close all;% Close all figures

path = 'E:\Data\BrainSlices\GluR2_antibody\181101_antiGluA2-CF633\3\4\';
outputpath='E:\Data\BrainSlices\GluR2_antibody\181101_antiGluA2-CF633\3\';
pathsplit=strsplit(path,'\');
outputindex=pathsplit(end-1);
outputname='_mem10_good50_maxparam1000'; % this will go to all the output filex e.g. result_tracking_picasso.txt
% frame selection
timerange=0; %if you want to specify frame = 1, otherwise = 0.
mint=1;
maxt=3000;
% ROI setting in nm.
useROI=0; % if you want to specify ROI 1, otherwise 0.
xmin=1632;
xmax=35520;
ymin=1344;
ymax=33216;
% parameters for z threshold and factor
threshold=600; 
factor=0.79; % for picasso if you have already multiplied the factor in localization step, this should be 1 for ThunderSTORM and qPALM  =0.79
% params for tracking analysis
param.mem=10;
param.dim=3; %for 3D use 3
param.good=50; %original 50
param.quiet=0;
maxdisp=1000;
dt=0.05; %exposure time. Default time step 0.05s=50ms
%for cluster analysis
k=50;
Eps=500;
save([path 'parameters_' outputname '.mat'], 'threshold','factor','param','maxdisp','dt','k','Eps');
%% Step1

if ~exist('fileName1','var')|| isempty(fileName1)
    [userfilein, userdirin]=uigetfile({
         '*.txt','Data file (*.txt)';'*.xlsx','Data file (*.xlsx)';...
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
         '*.txt','Data file (*.txt)';'*.xlsx','Data file (*.xlsx)';...
        '*.*','All Files (*.*)'},'Select the QD file to process',...
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

% if ~exist('threshold','var')||isempty(threshold)||~exist('factor','var')||isempty(factor)
%     prompt = {'Enter threshold:','Enter factor of z:'};
%     dlg_title = 'Input for preparing data';
%     num_lines = 1;
%     defaultanswer={'600','0.79'};
%     answer = inputdlg(prompt,dlg_title,num_lines,defaultanswer);
%     threshold=str2double(answer(1));
%     factor=str2double(answer(2));
% end

[path1,name1,ext1]=fileparts(fileName1);
[path2,name2,ext2]=fileparts(fileName2);
switch lower(ext1)
    case '.txt'
        data_PALM=textread(fileName1);
    case '.xlsx'
        data_PALM=xlsread(fileName1);
    otherwise
        error('Unexpected file extension for PALM: %s', ext1);
end

data_PALM(isnan(data_PALM))=Inf;
%%% Anne's revision
data_PALM=data_PALM(abs(data_PALM(:,7))<threshold,:);
data_PALM(:,7)=data_PALM(:,7)*factor;
%%%%%%%%%%%%%% revision
fid=fopen([userdirin 'data_PALM' outputname '.txt'],'w');
fprintf(fid,'%d %f %f %f %f %f %f\n',[data_PALM(:,1) data_PALM(:,2:7)]');    
fclose(fid);
%%%%%%%%%%%%%% revision
switch lower(ext2)
    case '.txt'
        data_QD=textread(fileName2);
    case '.xlsx'
        data_QD=xlsread(fileName2);
    otherwise
        error('Unexpected file extension for QD: %s', ext2);
end

data_QD(isnan(data_QD))=Inf;
data_QD=data_QD(abs(data_QD(:,7))<threshold,:);
data_QD(:,7)=data_QD(:,7)*factor;
data_QD=sortrows(data_QD,15);


if timerange ==1;
    idx_t=(data_QD(:,15)>=mint)&(dataQD(:,15)<=maxt);
    indices=find(idx);
    for i=(1:length(indices))
         dataQD_temp(i,:)=dataQD(indices(i),:);
    end
    dataQD=dataQD_temp;
else
end


%% Step2
x=data_QD(:,5);
y=data_QD(:,6);
z=data_QD(:,7);
t=data_QD(:,15);
positionlist=[x y z t];


result_tracking = track( positionlist, maxdisp, param );
%%%%%%%%%%%%%% revisionh
fid=fopen([userdirin 'result_tracking_' outputname '.txt'],'w');
fprintf(fid,'%f %f %f %d %d\n',result_tracking');    
fclose(fid);
%%%%%%%%%%%%%% revision
%ORIGINAL save('result_tracking','-ascii');
clear x y z t positionlist param maxdisp


%This m_file seperates particles according to their ID in the 'result_tracking', and
%calculates diffusion coeficient for each particle by fitting straight
%lines to the first 4 points of MSD/dt of the particle.
%of the tracking file and outputs the submatrices, e.g. sub_k is the
%submatrix of particle# k.
%load file
%make variable x y z frame number and particle ID
tic
x=result_tracking(:,1);
y=result_tracking(:,2);
z=result_tracking(:,3); % 0.79 is the correction factor for RI mismatch for water/oil for an 1.40 NA objective
frame=result_tracking(:,4);
particle=result_tracking(:,5);
d=length(particle); %total row number
ptotal=particle(d);%total particle number
tSteps=10;
clear Dif
clear MSD_Ind
figure
for j=1:ptotal
    dp=find(particle==j); %dimension of matrix of particle# j
    dplength=length(dp);
    dpmax=max(dp); % for particlej, the first row number
    dpmin=min(dp); % for particlej, the last row number
    m(j)=j; %test j value
    x1=x(dpmin:dpmax);
    y1=y(dpmin:dpmax);
    z1=z(dpmin:dpmax);
    frame1=frame(dpmin:dpmax);
    particle1=particle(dpmin:dpmax); % Added by Anne 1/14/2014
    TraceAll{j}=[x1,y1,z1,frame1,particle1]; % Modified by Anne 1/14/2014
    Trace_XYZ=[x1,y1,z1];
    % Anne added for range of trace
    TRI_tr = DelaunayTri(x1,y1,z1);
    [ch_tr v_tr] = convexHull(TRI_tr);
    TRI_tr2d = DelaunayTri(x1,y1);
    [ch_tr2d A_tr] = convexHull(TRI_tr2d);
    %[CtrT,SigT] = subclust(Trace_XYZ,1);
    Trace_range(j,1)=v_tr;
    Trace_range(j,2)=A_tr;
    %Trace_range(j,3)=SigT;
    %
    %eval(['sub_' num2str(j) '=[x1,y1,z1,frame1]']);% creating submatrix for particle# j sub_j
    [MSD00,d2r0,counts]=fMSD_vect(x1,y1,z1,frame1,dpmax,dpmin,tSteps); % Call function MSD, MSD.m file must be in the same folder as this file
    cutoff=10;  %cutoff for data points
    ind=find(counts>cutoff-1); % find index of counts that is equate to and above cutoff
    MSDCF=MSD00(ind); % find the MSD for those index from last line
    %ind1=[0,ind]; % Add (0,0) as the first point of the curve
    %MSDCF1=[0,MSDCF]; % Add (0,0) as the first point of the curve
    % ind1=ind(1:4);
    % MSDCF1=MSDCF(1:4);
    indlength=length(ind);
    if indlength>=4
        ind1=ind(1:4);
        MSDCF1=MSDCF(1:4);
        Dif(j,:)= polyfit(ind1',MSDCF1',1);
        
% ORIGINAL        f = fit(ind1',MSDCF1','poly1');
% ORIGINAL       Dif(j,:)=coeffvalues(f); % First column is the diffusion coeff (slope). To get um^2/s multiply by 2*10^(-5), to get cm^2/s, multiply by 2*10^(-13)
    else Dif(j,:)=[0,0]; ind1=ind; MSDCF1=MSDCF; %%% Changed on 2013/05/13, Dif has to have 2 columns, and if indlenth<4, ind1 and MSDCF1 has to have values
    end
    particle2(j)=j; % Added by Anne, 1/14/2014
    plot(ind1,MSDCF1)
    hold all
    %eval(['MSD00_' num2str(j) '=[MSD00]']);
    %eval(['d2r0_' num2str(j) '=[d2r0]']);
    %         clear dp dplength dpmax dpmin x1 y1 z1 frame1
    inst_msd{j}= inst_MSD( x1,y1,z1,frame1);
end
hold off
saveas(gcf,strcat(userdirin,['MSD_' outputname]),'fig')



Dif1=Dif(:,1);
Dif_positive=Dif1(Dif1>0);
Dif_time=Dif_positive/(dt*2*3); %D for time step dt, for 3D
Dif_track=Dif_positive/(dt*2*3*10^6);
Dif_track_ID=particle2(Dif1>0); % Added by Anne, 1/14/2014, for trace ID of corresponding D
Dif_track_W_ID(:,1)=Dif_track;
Dif_track_W_ID(:,2)=Dif_track_ID;
figure;
hDiff=histogram(log10(Dif_track),'NumBins',30,'Normalization','pdf');
xlabel('Diffusion coefficient (\mum^2/s, log)');
ylabel('Frequency');
%hist(log10(Dif_track), 30);
saveas(gcf,strcat(userdirin,['Dif_track_' outputname]),'fig')
save(strcat(userdirin,['Inst_diffusion_coefficients_' outputname '.mat']),'inst_msd');
save(strcat(userdirin,['diffusion_coefficients_' outputname '.mat']),'Dif_track_W_ID');
%YY commented save(strcat(userdirin,'Trace_range.mat'),'Trace_range')
clear x y z frame particle j d m x1 y1 z1 frame1 dp


%% Step3
%%%%%%%%%%%%%%%%%
%%%%%%%%Synapse_2
palmX=data_PALM(:,5);
palmY=data_PALM(:,6);
palmZ=data_PALM(:,7);
data_Syn=[palmX, palmY, palmZ];

% if ~exist('k','var')||isempty(k)||~exist('Eps','var')||isempty(Eps)
%     prompt = {'Enter k - number of nearby points:','Enter Eps - neighborhood range:'};
%     dlg_title = 'Input for synapse analysis';
%     num_lines = 1;
%     defaultanswer={'50','500'};
%     answer = inputdlg(prompt,dlg_title,num_lines,defaultanswer);
%     k=str2double(answer(1));
%     Eps=str2double(answer(2));
% end

% Finding clusters
[Class,type]=dbscan_conservative(data_Syn,k,Eps); 
% Make new matrix  
Syn=[palmX,palmY,palmZ,Class',type'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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



%% trace_range.m This script calculates the range of traces in x, y and z direction
clearvars -except path outputname outputpath outputindex;%clear
%close all;% Close all figures

% %% Read file
% if ~exist('fileName1','var')|| isempty(fileName1)
%     [userfilein, userdirin]=uigetfile({
%         '*.txt','Data file (*.txt)';...
%         '*.*','All Files (*.*)'},'Select the result_tracking file to process',...
%         'E:\Data\BidentateQD\NewOHsQD\190607_SA-Atto647N_monothiolQD_temperature\3\4\');
%     fileName1=fullfile(userdirin,userfilein);
% else
%     if ~exist(fileName1,'file')
%         fprintf('File not found: %s\n',fileName1);
%         return;
%     end
% end
fileName1=[path 'result_tracking_' outputname '.txt'];
result=textread(fileName1); % Tracking result



%% Group trace according to ID and arrange in cells (unit nm)
x=result(:,1);
y=result(:,2);
z=result(:,3);
frame=result(:,4);
particle=result(:,5);
d=length(particle); %total row number
ptotal=particle(d);%total particle number
for k=1:ptotal
    clear x1 y1 z1 frame1 particle1 dp 
        dp=find(particle==k); %dimension of matrix of particle# K 
        dplength=length(dp);
        dpmax=max(dp);
        dpmin=min(dp);
        m(k)=k; %test k value
        x1=x(dpmin:dpmax);
        y1=y(dpmin:dpmax);
        z1=z(dpmin:dpmax);
        frame1=frame(dpmin:dpmax);
        particle1=particle(dpmin:dpmax);
        trace{k}=[x1,y1,z1,frame1,particle1]; 
        trace_R(k,1)=max(x1)-min(x1); % trace range in X
        trace_R(k,2)=max(y1)-min(y1); % trace range in Y
        trace_R(k,3)=max(z1)-min(z1); % trace range in z
        trace_R(k,4)=max(frame1)-min(frame1);
        trace_R(k,5)=k;    % trace ID       
end
figure
trace_R(:,6)=max(trace_R(:,1:3)');
hTrajectory=histogram(log10(trace_R(:,6)/1000),48,'Normalization','pdf');
%hist(log10(trace_R(:,6)/1000),48) % trace_r in um
xlabel('Trajectory length (\mum, log)');
ylabel('Frequency');
userdirin=path;
save(strcat(userdirin,['trace_R_' outputname '.mat']),'trace_R'); % trace_R: [X Y Z frame ID max(X Y Z)]
saveas(gcf,strcat(userdirin,['hist_traceRange_' outputname]),'fig');
clear x y z frame particle k d m x1 y1 z1 frame1 dp ptotal dplength dpmax dpmin

%% trace_diffusion.m
clearvars -except path outputname outputpath outputindex;%clear
% close all;% Close all figures
% 


% if ~exist('fileName1','var')|| isempty(fileName1)
%     [userfilein, userdirin]=uigetfile({
%         '*.mat','Data file (*.mat)';...
%         '*.*','All Files (*.*)'},'Select the trace_R file to process',...
%         path);
%     fileName1=fullfile(userdirin,userfilein);
% else
%     if ~exist(fileName1,'file')
%         fprintf('File not found: %s\n',fileName1);
%         return;
%     end
% end
fileName1=[path 'trace_R_' outputname '.mat'];
% if ~exist('fileName2','var')|| isempty(fileName2)
%     [userfilein, userdirin]=uigetfile({
%         '*.mat','Data file (*.mat)';...
%         '*.*','All Files (*.*)'},'Select the diffusion_coefficient file to process',...
%         path);
%     fileName2=fullfile(userdirin,userfilein);
% else
%     if ~exist(fileName2,'file')
%         fprintf('File not found: %s\n',fileName2);
%         return;
%     end
% end
fileName2=[path 'diffusion_coefficients_' outputname '.mat'];

% path = 'E:\experiment\20161202_NMDAR_TRACKING\SA\1nM\FM\1\fPALM\diffusion\';
% 
% load('E:\experiment\20161202_NMDAR_TRACKING\SA\1nM\FM\1\fPALM\diffusion\trace_R.mat')
% load('E:\experiment\20161202_NMDAR_TRACKING\SA\1nM\FM\1\fPALM\diffusion\diffusion_coefficients.mat')
% tr = trace_R;
% s = size(tr);
% s1 = s(1);
% 
% d = Dif_track_W_ID;
% p = size(d);
% p1 = p(1);

% path = 'F:\experiment\NR2A+NR2B_tracking\20180816_DIV15_NR2A_SA_Atto647N\1\test between perisynaptic and synaptic\';
% file1 = 'trace_R.mat';
% file2 = 'diffusion_coefficients_605.mat';
% fname1 = fullfile(path,file1);
% fname2 = fullfile(path,file2);
load(fileName1);
tr = trace_R;
s = size(tr);
s1 = s(1);

load(fileName2);
d = Dif_track_W_ID;
p = size(d);
p1 = p(1);
%
output= zeros(p1,2);
output2= zeros(p1,4);

for k=1:p1;
   
    output(k,:) = [tr(d(k,2),6) d(k,1)];

    
end


output2(:,1)=output(:,1);
output2(:,2)=output(:,2);

output2(:,3) = log10(output(:,1)/1000);
output2(:,4) = log10(output(:,2));
userfilein=['diffusion_coefficients' outputname '.mat'];
userdirin=path;
userfilein_split  = strsplit(userfilein, '.');
matched = strcat(userdirin, userfilein_split(1), '-diffusion_trace.xlsx');
if exist('output2') ==1
    xlswrite(matched{1}, output2);
    xlswrite([outputpath userfilein_split{1} '-diffusion_trace_' outputindex{1} '.xlsx'], output2);
else
end

%% TraceDisperse_modified.m
clearvars -except path outputname outputpath outputindex;%clear
% close all;% Close all figures

% if ~exist('fileName1','var')|| isempty(fileName1)
%     [userfilein, userdirin]=uigetfile({
%         '*.txt','Data file (*.txt)';...
%         '*.*','All Files (*.*)'},'Select the result_tracking to process',...
%         path);
%     fileName1=fullfile(userdirin,userfilein);
% else
%     if ~exist(fileName1,'file')
%         fprintf('File not found: %s\n',fileName1);
%         return;
%     end
% end
fileName1=[path 'result_tracking_' outputname '.txt'];
% if ~exist('fileName2','var')|| isempty(fileName2)
%     [userfilein, userdirin]=uigetfile({
%         '*.mat','Data file (*.mat)';...
%         '*.*','All Files (*.*)'},'Select the Syanpse1 file to process',...
%         path);
%     fileName2=fullfile(userdirin,userfilein);
% else
%     if ~exist(fileName2,'file')
%         fprintf('File not found: %s\n',fileName2);
%         return;
%     else [userdirin,~,~]=fileparts(fileName2);
%         userdirin=strcat(userdirin,'\');
%     end
% end
fileName2=[path 'Synapses1_' outputname '.mat'];
result_tracking=textread(fileName1);
S_Synapses1=load(fileName2);
Synapses1=S_Synapses1.Synapses1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step2
%%%%%%%%%%%%%
%%%%%%trace_center
% This script finds the center of traces abtained from tracking.m
%===============================
% load data from tracking.m
x=result_tracking(:,1);
y=result_tracking(:,2);
z=result_tracking(:,3);
frame=result_tracking(:,4);
particle=result_tracking(:,5);
d=length(particle); %total row number
ptotal=particle(d); %total particle number
for n=1:ptotal
    dp=find(particle==n); %dimension of matrix of particle# K
%     dplength=length(dp);
    dpmax=max(dp);
    dpmin=min(dp);
    m(n)=n; %test k value
    x1=x(dpmin:dpmax);
    y1=y(dpmin:dpmax);
    z1=z(dpmin:dpmax);
    pt1=particle(dpmin:dpmax);
    %trace=[x1,y1,z1];
    frame1=frame(dpmin:dpmax);
    traces{n}=[x1,y1,z1,frame1,pt1];% creating a cell array for particle# k traces{k}
    traces_center(n,1)=mean(x1);
    traces_center(n,2)=mean(y1);
    traces_center(n,3)=mean(z1);
    

% save(strcat(userdirin,'traces.mat'),'traces');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  nSyn=length(Synapses1);
  
  % syn_dis=2000;  % unit in nm
  
 
    
      for rg=1:nSyn
        syn_rg=Synapses1{rg}; % load data
        syn_length=length(syn_rg);
        % Find center of the synapse
        clear trace_nearby syn_traces
        %radii=1; % Edit by Anne 11/26/2013 
        Ctr0=mean(syn_rg);%ORIGINAL         [Ctr00,Sig00]=subclust(syn_rg,radii); % Find center of a subcluster and the range of influence of the center
        %%make sure Ctr0 has only 1 row
        dist=zeros(length(traces{1}));
        % Find traces close to the synapse
        syn_traces(rg,1)=(Ctr0(1)-traces_center(n,1)).*(Ctr0(1)-traces_center(n,1));
        syn_traces(rg,2)=(Ctr0(2)-traces_center(n,2)).*(Ctr0(2)-traces_center(n,2));
        syn_traces(rg,3)=(Ctr0(3)-traces_center(n,3)).*(Ctr0(3)-traces_center(n,3));
        distance(rg)=sqrt(syn_traces(rg,1)+syn_traces(rg,2)+syn_traces(rg,3))/1000;
        dist = [distance];
        dist2 = dist';
        mini=min(dist2);
        t= find(dist2 == mini);
        a=Synapses1{t};
        syn_cent0= mean(a);
        Syn_ctrX = syn_cent0(1);
        Syn_ctrY = syn_cent0(2);
        Syn_ctrZ = syn_cent0(3);
       
       
        TS_X=(traces{n}(:,1)-Syn_ctrX)/1000; % Trace to synapse distance delta X
        TS_Y=(traces{n}(:,2)-Syn_ctrY)/1000; 
        TS_Z=(traces{n}(:,3)-Syn_ctrZ)/1000;
          
        TS_dist{n}=sqrt(TS_X.*TS_X+TS_Y.*TS_Y+TS_Z.*TS_Z);
        %syn_nearby{rg}=find(syn_traces(:,4)<syn_dis);
        %syn_nearby_t{rg}=find(syn_traces(:,4)<syn_dis)'; % syn_nearby_t is for coverting into traceID array
        %Syn_ctr(rg,:)=Ctr0(1:3);% Add by Anne 11/26/2013
      end
    
   % again synapess center

      
end
userdirin=path;
 save(strcat(userdirin,['TS_dist_' outputname '.mat']),'TS_dist');
 DistBin=(0:0.02:2);
 figure;
 subplot(1,2,1);
 hDistance=histogram(cell2mat(TS_dist'),DistBin,'Normalization','pdf');
 %hist(cell2mat(TS_dist'),DistBin)
 xlim([0 1.98]);
 xlabel('Distance (\mum), PDF');
 ylabel('Frequency');
 subplot(1,2,2);
 hDistanceCDF=histogram(cell2mat(TS_dist'),DistBin,'Normalization','cdf');
 xlim([0 1.98]);
 xlabel('Distance (\mum), CDF');
 ylabel('CDF');
saveas(gcf,strcat(userdirin,['TS_dist_' outputname]),'fig');

distance=TS_dist';

D=cell2mat(distance);

% L=size(D_table);



% fclose(fid);
%xlswrite(strcat(path,['TS_dist_cal_' [outputname '_' outputindex{1}] '.xlsx']), D);
xlswrite(strcat(outputpath,['TS_dist_cal_' [outputname '_' outputindex{1}] '.xlsx']), D);
%     save(strcat(userdirin,'Syn_ctr.mat'),'Syn_ctr');
%     syn_nearby0=find(~cellfun(@isempty,syn_nearby)); % Synapses ID with traces around
  
  

