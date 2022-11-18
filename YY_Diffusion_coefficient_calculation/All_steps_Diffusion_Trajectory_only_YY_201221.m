%line 39 threshold fixed 600 nm 0.79
% line 181 enter k and eps: fixed 50, 500
% Modified by Anne. Comapred to V2, V3 used tracking fucntion on mTraceIN
% and mTraceOUT in order to avoid connecting discontinued tracdes of same particle in or
% out of synapse
% Modified by Anne, 1/14/2014, added Dif_track ID, for note which trace ID is for corrisponding Dif_track
%% Step 0 Parameters. YOU NEED TO TAKE A LOOK AT THIS PART CAREFULLY, OTHERWISE YOU WILL GET ERRORS OR MISINTERPRETED RESULTS.
clc;clear;
%close all;% Close all figures

path = 'H:\Data_from07222021\DNA-PAINT(Live-cell)\SPT_tracking\cLTP_tracking\211203_cLTPCTRL_DNA-PAINT_DIV16\';
outputpath=[path 'Results\'];
pathsplit=strsplit(path,'\');
outputindex=pathsplit(end-1);
outputname='postvehicle_after'; % this will go to all the output filex e.g. result_tracking_picasso.txt
isitpicasso=0;
actualdimension=2;
% frame selection
timerange=0; %if you want to specify frame = 1, otherwise = 0.
mint=1;
maxt=3000;
% ROI setting in nm.
useROI=0; % if you want to specify ROI 1, otherwise 0.
xmin=11236;
xmax=27326.8;
ymin=3392;
ymax=14755.2;
% parameters for z threshold and factor
threshold=500; %YY corrected for ThunderSTORM. originally 600 from line ~59
if isitpicasso==1;
    factor=1;
else
    factor=0.79; % for picasso if you have already multiplied the factor in localization step, this should be 1 for ThunderSTORM =0.79
end
% params for tracking analysis
param.mem=3;
param.dim=3; %for 3D use 3
param.good=5; %original 50
param.quiet=0;
maxdisp=500;
dt=0.1; %exposure time. Default time step 0.05s=50ms
tSteps=10; %length of MSD
diffusion_thres = 4; % how many time points used in calculating diffusion coefficients
cutoff=param.good; % in MSD calculation minimum length of trajectory line 249
%for cluster analysis
% k - number of objects in a neighborhood of an object
% (minimal number of objects considered as a cluster)
% Eps - neighborhood radius, if not known avoid this parameter or put []
filter_photon=0; % if you filter non-specific binding using photon number 1, otherwise 0.
photon_cutoff=6000;

%% Step1

if ~exist('fileName2','var')|| isempty(fileName2)
    [userfilein, userdirin]=uigetfile({
         '*.csv','Data file (*.csv)';...
        '*.*','All Files (*.*)'},'Select the AMPAR(or tracking) file to process',...
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
% threshold=1000; %YY corrected for ThunderSTORM. originally 600
% factor=0.79;

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


%% Step2 tracking.
if actualdimension==3;
    x=data_QD(:,3);
    y=data_QD(:,4);
    z=data_QD(:,5);% z=data_QD(:,5); %for 3D
    t=data_QD(:,2);
    positionlist=[x y z t]; %for 3D use this [x y z t];
elseif actualdimension==2;
    x=data_QD(:,3);
    y=data_QD(:,4);
    z=zeros(length(x),1);
    t=data_QD(:,2);
    positionlist=[x y z t];
end



%result_tracking = track( positionlist, maxdisp, param );
if actualdimension==3;
    [result_tracking photon_traj] = track_YY(positionlist,data_QD(:,8),maxdisp,param);
elseif actualdimension ==2;
    [result_tracking photon_traj] = track_YY(positionlist,data_QD(:,6),maxdisp,param);
end
    
ID=[];
if filter_photon==1;
    for i=1:max(result_tracking(:,5));
        if any(photon_traj(find(result_tracking(:,5)==i)) > photon_cutoff);
            ID=[ID;i];
        end
    end
    filtered_tracking=[];
    for i=1:length(ID);
        ID_temp=find(result_tracking(:,5)==ID(i,1));
        filtered_tracking=[filtered_tracking; result_tracking(ID_temp,:)];
    end
    result_tracking=filtered_tracking;
end
% need to rearrange ID       
acc=unique(result_tracking(:,5));
for i=1:length(result_tracking);
    result_tracking(i,5)=find(result_tracking(i,5)==acc);
end

%%%%%%%%%%%%%% revisionh
fid=fopen([outputpath 'result_tracking_' outputname '.txt'],'w');
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
    z=result_tracking(:,3); % for 3D 0.79 is the correction factor for RI mismatch for water/oil for an 1.40 NA objective
    frame=result_tracking(:,4); 
    particle=result_tracking(:,5);
    d=length(particle); %total row number
    ptotal=particle(d);%total particle number
    %tSteps=10;
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
    % YY means Yeoan's edition
    % Anne added for range of trace
    %YY TRI_tr = DelaunayTri(x1,y1,z1);
    %YY [ch_tr v_tr] = convexHull(TRI_tr);
    %YY TRI_tr2d = DelaunayTri(x1,y1);
    %YY [ch_tr2d A_tr] = convexHull(TRI_tr2d);
    %[CtrT,SigT] = subclust(Trace_XYZ,1);
    %YY Trace_range(j,1)=v_tr;
    %YY Trace_range(j,2)=A_tr;
    %Trace_range(j,3)=SigT;
    %
    %eval(['sub_' num2str(j) '=[x1,y1,z1,frame1]']);% creating submatrix for particle# j sub_j
        [MSD00,d2r0,counts]=fMSD_vect(x1,y1,z1,frame1,dpmax,dpmin,tSteps); % Call function MSD, MSD.m file must be in the same folder as this file
        %cutoff=10;  %cutoff for data points YY
        ind=find(counts>cutoff-1); % find index of counts that is equate to and above cutoff
        MSDCF=MSD00(ind); % find the MSD for those index from last line
    %ind1=[0,ind]; % Add (0,0) as the first point of the curve
    %MSDCF1=[0,MSDCF]; % Add (0,0) as the first point of the curve
    % ind1=ind(1:4);
    % MSDCF1=MSDCF(1:4);
        indlength=length(ind);
        if indlength>=diffusion_thres;
            ind1=ind(1:diffusion_thres);
            MSDCF1=MSDCF(1:diffusion_thres);
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
    saveas(gcf,strcat(outputpath,['MSD_' outputname]),'fig')

Dif1=Dif(:,1);
Dif_positive=Dif1(Dif1>0);
if actualdimension==3;
    Dif_time=Dif_positive/(dt*2*3); %D for time step dt, for 3D
    Dif_track=Dif_positive/(dt*2*3*10^6);
else actualdimension==2;
        Dif_time=Dif_positive/(dt*2*2); %D for time step dt, for 2D
    Dif_track=Dif_positive/(dt*2*2*10^6);
end
Dif_track_ID=particle2(Dif1>0); % Added by Anne, 1/14/2014, for trace ID of corresponding D
Dif_track_W_ID(:,1)=Dif_track;
Dif_track_W_ID(:,2)=Dif_track_ID;

figure;
hDiff=histogram(log10(Dif_track),'NumBins',30,'Normalization','pdf');
xlim([-5 0]);
xlabel('Diffusion coefficient (\mum^2/s, log)');
ylabel('Frequency');
%hist(log10(Dif_track), 30);
saveas(gcf,strcat(outputpath,['Dif_track_' outputname]),'fig')
save(strcat(outputpath,['Inst_diffusion_coefficients_' outputname '.mat']),'inst_msd');
save(strcat(outputpath,['diffusion_coefficients_' outputname '.mat']),'Dif_track_W_ID');
%YY commented save(strcat(userdirin,'Trace_range.mat'),'Trace_range')
clear x y z frame particle j d m x1 y1 z1 frame1 dp


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
fileName1=[outputpath 'result_tracking_' outputname '.txt'];
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
        for j=1:length(x1);
            distances=sqrt((x1-x1(j)).^2+(y1-y1(j)).^2+(z1-z1(j)).^2);
        end
        trace_R(k,1)=max(x1)-min(x1); % trace range in X
        trace_R(k,2)=max(y1)-min(y1); % trace range in Y
        trace_R(k,3)=max(z1)-min(z1); % trace range in z
        trace_R(k,4)=max(frame1)-min(frame1);
        trace_R(k,5)=k;    % trace ID
        trace_R(k,6)=max(distances);
end
figure
% trace_R(:,6)=max(trace_R(:,1:3)');
hTrajectory=histogram(log10(trace_R(:,6)/1000),48,'Normalization','pdf');
%hist(log10(trace_R(:,6)/1000),48) % trace_r in um
xlabel('Trajectory length (\mum, log)');
ylabel('Frequency');
userdirin=path;
save(strcat(outputpath,['trace_R_' outputname '.mat']),'trace_R'); % trace_R: [X Y Z frame ID max(X Y Z)]
saveas(gcf,strcat(outputpath,['hist_traceRange_' outputname]),'fig');
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
fileName1=[outputpath 'trace_R_' outputname '.mat'];
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
fileName2=[outputpath 'diffusion_coefficients_' outputname '.mat'];

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
    %xlswrite(matched{1}, output2);
    xlswrite([outputpath userfilein_split{1} '-diffusion_trace_' outputindex{1} '.xlsx'], output2);
else
end



