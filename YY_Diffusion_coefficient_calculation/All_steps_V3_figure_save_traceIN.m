function All_steps(fileName1,fileName2,threshold,factor,k,Eps,optional)
% Modified by Anne. Comapred to V2, V3 used tracking fucntion on mTraceIN
% and mTraceOUT in order to avoid connecting discontinued tracdes of same particle in or
% out of synapse
%% Step1
clc;%clear
close all;% Close all figures
if ~exist('fileName1','var')|| isempty(fileName1)
    [userfilein, userdirin]=uigetfile({
        '*.xls','Data file (*.xls)';...
        '*.*','All Files (*.*)'},'Select the PALM file to process',...
        'D:\mydocument\MATLAB\PALM_tracking_AnalysisCode');
    fileName1=fullfile(userdirin,userfilein);
else
    if ~exist(fileName1,'file')
        fprintf('File not found: %s\n',fileName1);
        return;
    end
end
if ~exist('fileName2','var')|| isempty(fileName2)
    [userfilein, userdirin]=uigetfile({
        '*.xls','Data file (*.xls)';...
        '*.*','All Files (*.*)'},'Select the QD file to process',...
        'D:\mydocument\MATLAB\PALM_tracking_AnalysisCode');
    fileName2=fullfile(userdirin,userfilein);
else
    if ~exist(fileName2,'file')
        fprintf('File not found: %s\n',fileName2);
        return;
    else [userdirin,~,~]=fileparts(fileName2);
        userdirin=strcat(userdirin,'\');
    end
end

if ~exist('threshold','var')||isempty(threshold)||~exist('factor','var')||isempty(factor)
    prompt = {'Enter threshold:','Enter factor of z:'};
    dlg_title = 'Input for preparing data';
    num_lines = 1;
    defaultanswer={'600','0.79'};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultanswer);
    threshold=str2double(answer(1));
    factor=str2double(answer(2));
end
data_PALM=xlsread(fileName1);
data_PALM(isnan(data_PALM))=Inf;
%%% Anne's revision
data_PALM=data_PALM(abs(data_PALM(:,7))<threshold,:);
data_PALM(:,7)=data_PALM(:,7)*factor;
%%%%%%%%%%%%%% revision
fid=fopen([userdirin 'data_PALM.txt'],'w');
fprintf(fid,'%d %f %f %f %f %f %f\n',[data_PALM(:,1) data_PALM(:,2:7)]');    
fclose(fid);
%%%%%%%%%%%%%% revision
data_QD=xlsread(fileName2);
data_QD(isnan(data_QD))=Inf;
data_QD=data_QD(abs(data_QD(:,7))<threshold,:);
data_QD(:,7)=data_QD(:,7)*factor;
data_QD=sortrows(data_QD,15);
%% Step2
x=data_QD(:,5);
y=data_QD(:,6);
z=data_QD(:,7);
t=data_QD(:,15);
positionlist=[x y z t];
param.mem=10;
param.dim=3;
param.good=30;
param.quiet=0;
maxdisp=500;
result_tracking = track( positionlist, maxdisp, param );
%%%%%%%%%%%%%% revision
fid=fopen([userdirin 'result_tracking.txt'],'w');
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
    TraceAll{j}=[x1,y1,z1,frame1];
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
    plot(ind1,MSDCF1)
    hold all
    %eval(['MSD00_' num2str(j) '=[MSD00]']);
    %eval(['d2r0_' num2str(j) '=[d2r0]']);
    %         clear dp dplength dpmax dpmin x1 y1 z1 frame1
    inst_msd{j}= inst_MSD( x1,y1,z1,frame1);
end
hold off
dt=0.05; %Default time step 0.05s=50ms
Dif1=Dif(:,1);
Dif_positive=Dif1(Dif1>0);
Dif_time=Dif_positive/(dt*2*3); %D for time step dt, for 3D
Dif_track=Dif_positive/(dt*2*3*10^6);
figure
hist(Dif_track, 30);
saveas(gcf,strcat(userdirin,'Dif_track'),'fig')
save(strcat(userdirin,'Inst_diffusion_coefficients.mat'),'inst_msd');
save(strcat(userdirin,'diffusion_coefficients.mat'),'Dif_track');
save(strcat(userdirin,'Trace_range.mat'),'Trace_range')
clear x y z frame particle j d m x1 y1 z1 frame1 dp

%% Step3
%%%%%%%%%%%%%%%%%
%%%%%%%%Synapse_2
palmX=data_PALM(:,5);
palmY=data_PALM(:,6);
palmZ=data_PALM(:,7);
data_Syn=[palmX, palmY, palmZ];

if ~exist('k','var')||isempty(k)||~exist('Eps','var')||isempty(Eps)
    prompt = {'Enter k - number of nearby points:','Enter Eps - neighborhood range:'};
    dlg_title = 'Input for synapse analysis';
    num_lines = 1;
    defaultanswer={'50','500'};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultanswer);
    k=str2double(answer(1));
    Eps=str2double(answer(2));
end

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
saveas(gcf,strcat(userdirin,'syn'),'fig');
save(strcat(userdirin,'Synapses1.mat'),'Synapses1');
axis equal


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%Syn_size
% This script calculates the size of synapses found by synapse_2.m. X, Y ,Z
% and volumn is calculated
[synNumb,synNull]=size(Synapses1); % get number of synapses
% SynapseSize=zeros(synNumb,5);
Numb=[1:synNumb];
sizeX=zeros(synNumb,1);
sizeY=zeros(synNumb,1);
sizeZ=zeros(synNumb,1);
Vol=zeros(synNumb,1);
Number=zeros(synNumb,1);
Density=zeros(synNumb,1);
figure
%Volume and size
for g=1:synNumb
    m=Synapses1{g};
    SynX=m(:,1);
    SynY=m(:,2);
    SynZ=m(:,3);
    %Volume
    TRI = DelaunayTri(SynX,SynY,SynZ);
    [ch v] = convexHull(TRI);
    Vol(g)=v;
    K{g}=ch;
    hold all
    trisurf(ch, TRI.X(:,1),TRI.X(:,2),TRI.X(:,3))
    %Area
    %     DT_XY = DelaunayTri(SynX,SynY);
    %     [cha1 a1] = convexHull(DT_XY);
    %     Area_XY(g)=a1;
    %     DT_YZ = DelaunayTri(SynY,SynZ);
    %     [cha2 a2] = convexHull(DT_YZ);
    %     Area_YZ(g)=a2;
    %     DT_XZ = DelaunayTri(SynX,SynZ);
    %     [cha3 a3] = convexHull(DT_XZ);
    %     Area_XZ(g)=a3;
    %linear size
    sizeX(g)=max(SynX)-min(SynX);
    sizeY(g)=max(SynY)-min(SynY);
    sizeZ(g)=max(SynZ)-min(SynZ);
    Number(g)=size(Synapses1{g},1);
    Density(g)=Number(g)/Vol(g);
end
SynapseSize=[Numb', Vol,sizeX, sizeY, sizeZ,Number,Density];
saveas(gcf,strcat(userdirin,'syn_D'),'fig');
save(strcat(userdirin,'SynapseSize.mat'),'SynapseSize');
hold off

%% Step4
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
end
save(strcat(userdirin,'traces.mat'),'traces');
clear x y z frame particle k d m x1 y1 z1 frame1 dp ptotal dpmax dplength dpmin pt1


%%%%%%%%%%%%%
%%%%%%%%nearby_trace
if ~exist('optional','var')|| isempty(optional)
    optionalAnswer=questdlg('Do you want to run the Nearby_trace?','Save Mat?','yes','no','yes');
    if strcmp('yes',optionalAnswer)
        optional=1;
    else
        optional=0;
    end
end
% This script detects nearby traces based on the center of the current
% synapse.
% Need to load Synapses and run trace_center.m before this script.
clear syn_traces trace_n_syn1 syn_w_trace1 k kl mk mk_l Ctr0 Ctr00 Sig0 Sig00 nSyn nSyn2 radii rg syn_dis syn_length syn_rg syn_sub
if optional
    %===================================
    % Load synapse
    %load('Synapses.mat')
    nSyn=length(Synapses1);
    % Load trace and find the centers of the traces
    %load('traces_center.mat')
    %===========================
    % Define a range around synapse center to find traces that is close to the
    % synapse
    syn_dis=2000; % unit in nm
    %===========================
    % Loop
    for rg=1:nSyn
        syn_rg=Synapses1{rg}; % load data
        syn_length=length(syn_rg);
        % Find center of the synapse
        clear trace_nearby syn_traces
        %radii=1; % Edit by Anne 11/26/2013 
        Ctr0=mean(syn_rg);%ORIGINAL         [Ctr00,Sig00]=subclust(syn_rg,radii); % Find center of a subcluster and the range of influence of the center
        %%make sure Ctr0 has only 1 row
        %         if length(Ctr00(:,1))==1
        %             Ctr0=Ctr00;
        %         else
        %             Ctr0=mean(Ctr00);
        %         end
        %         if length(Sig00(:,1)==1)
        %             Sig0=Sig00;
        %         else
        %             Sig0=mean(Sig00);
        %         end
        %         syn_sub(rg,:)=Ctr0(1:3);
        % Find traces close to the synapse
        syn_traces(:,1)=(traces_center(:,1)-Ctr0(1)).*(traces_center(:,1)-Ctr0(1));
        syn_traces(:,2)=(traces_center(:,2)-Ctr0(2)).*(traces_center(:,2)-Ctr0(2));
        syn_traces(:,3)=(traces_center(:,3)-Ctr0(3)).*(traces_center(:,3)-Ctr0(3));
        syn_traces(:,4)=sqrt(syn_traces(:,1)+syn_traces(:,2)+syn_traces(:,3));
        
        syn_nearby{rg}=find(syn_traces(:,4)<syn_dis);
        %syn_nearby_t{rg}=find(syn_traces(:,4)<syn_dis)'; % syn_nearby_t is for coverting into traceID array
        Syn_ctr(rg,:)=Ctr0(1:3);% Add by Anne 11/26/2013
    end
    save(strcat(userdirin,'Syn_ctr.mat'),'Syn_ctr');
    syn_nearby0=find(~cellfun(@isempty,syn_nearby)); % Synapses ID with traces around
    
    % Plot
    nSyn2=length(syn_nearby0);
    figure
    hold all
    % for rg=1:nSyn
    %     syn_rg=Synapses{rg};
    %     scatter3(syn_rg(:,1), syn_rg(:,2), syn_rg(:,3),30,'filled','c')
    % end
    %%Plot all tracse in the figure
    clear trace_K traceID
    traceID_length=size(traces');
    for traceID=1:traceID_length
        trace_K=traces{traceID};
        plot3(trace_K(:,1),trace_K(:,2),trace_K(:,3),'g')
    end

    %clear trace_K traceID
   clear Syn_nearby_AveR
    for k=1:nSyn2
        clear mk mk_l
        mk=syn_nearby{syn_nearby0(k)}; % syn_nearby
        mk_l=length(mk);
        for kl=1:mk_l
            clear mk_trace mk_trace1 % add by Anne on 11/26/2013
            mk_trace=traces{mk(kl)};
            plot3(mk_trace(:,1),mk_trace(:,2),mk_trace(:,3),'r');
            
            % Add by Anne on 12/15/2013
            Syn_ctrX=Syn_ctr(syn_nearby0(k),1);
            Syn_ctrY=Syn_ctr(syn_nearby0(k),2);
            Syn_ctrZ=Syn_ctr(syn_nearby0(k),3);
            TS_X=mk_trace(:,1)-Syn_ctrX; % Trace to synapse distance delta X
            TS_Y=mk_trace(:,2)-Syn_ctrY; 
            TS_Z=mk_trace(:,3)-Syn_ctrZ;
            Ave_TS_dist(kl)=mean(sqrt(TS_X.*TS_X+TS_Y.*TS_Y+TS_Z.*TS_Z)/1000); % average distance from trace to synapse center 
            
            mk_trace1{kl}=mk_trace; % Add by Anne on 11/26/2013
        end
        Syn_nearby_AveR{syn_nearby0(k)}=Ave_TS_dist;
        clear Ave_TS_dist
        mk_syn=Synapses1{syn_nearby0(k)};
        trace_n_syn1{syn_nearby0(k)}=cell2mat(mk_trace1');
        scatter3(mk_syn(:,1), mk_syn(:,2), mk_syn(:,3),10,'filled','b');
        %%Output cell array of trace that is close to synapses (the following code is not correct)
        % the following code is not correct 
        %mk_trace1=mk_trace;
        %mk_trace2(:,6)=syn_nearby0(k)*ones(length(mk_trace),1);
        % for each cell, x,y,z,frame,particle% Edit by Anne 11/26/2013
        %%Output cell array of synapses with traces around
        %syn_w_trace1{k}=mk_syn; % for each cell, x,y,z,class,synapse
        %%Clear
        %clear mk_1 mk_trace mk_syn mk_trace1
        
        %Add by Anne 11/26/2013
        
         % Distance of each point on trace to the center of synapse in um
    end
    
    save(strcat(userdirin,'trace_n_syn1.mat'),'trace_n_syn1');% Add by Anne 11/26/2013
    save(strcat(userdirin,'Syn_nearby_AveR.mat'),'Syn_nearby_AveR');% Add by Anne 11/26/2013
    %trace_n_syn=trace_n_syn1';
    %syn_w_trace=syn_w_trace1';
    mTraceID=cell2mat(syn_nearby');
    TraceID=unique(mTraceID);
    
    hold off
    axis equal
    saveas(gcf,strcat(userdirin,'syn_trace'),'fig');
    clear syn_traces syn_w_trace1 k kl mk mk_l Ctr0 Ctr00 Sig0 Sig00 nSyn nSyn2 radii rg syn_dis syn_length syn_rg syn_sub trace_nearby
end
Ave_TS_dist=cell2mat(Syn_nearby_AveR); % add by Anne 12/15/2013
figure; hist(Ave_TS_dist,28) % add by Anne 12/15/2013
saveas(gcf,strcat(userdirin,'Ave_TS_dist'),'fig');
length(find(Ave_TS_dist<1))
length(Ave_TS_dist)
In_syn_trace_percent_1um=length(find(Ave_TS_dist<1))/length(Ave_TS_dist)
In_syn_trace_percent_750nm=length(find(Ave_TS_dist<0.75))/length(Ave_TS_dist)
In_syn_trace_percent_500nm=length(find(Ave_TS_dist<0.5))/length(Ave_TS_dist)


%% Step 5
if ~optional
    for n=1:size(Synapses1,1)
    syn_nearby{n}=(1:size(traces,2))';
    end
end
% Syn_COAR_all: For all the synapses, syanpse center, orientation, axis, range
% This script rotates a 3D cluster data (synapse) to its principal axis and
% calculates its center and sigma. A synaptic range is decided based on the
% ellipsoid fitted to the original scattered data. Than the synapse data
% and ellipsoid is rotated back to the original coordinate system.

% Get data
% load('Synapses.mat')
% syn_n=length(Synapses1);
% %rg=1;
% figure(6)
% hold all
% for rg=1:syn_n % Define which synapse
%     %%load data
%     syn_rg=Synapses1{rg}; % load data
%     if isempty(syn_nearby{rg})==0 % find traces near this synapse
%         syn_trace_rg=syn_nearby{rg}; % syn_trace_rg is the ID of traces that is near synapse rg.
%     else
%         syn_trace_rg=0;
%     end
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Principal analysis
%     % Find the principal axis. Keep the scattered data centeted, so that the
%     % new scattered data from PCA can be rotated easily back to the original
%     % coordinate system.
%     [SynR_coeff, SynR_score, latent]=pca(syn_rg(:,1:3),'Centered',false);
%     
%     syn_rr=SynR_score*inv(SynR_coeff); % syn_rr is the data rotated back to the original coordinate system.
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%Syn_angle
%     %   This script calculates the angle alpha between z axis and the major
%     %   axis of synapses
%     %%Calculate the angle between z-axis and the major axis of the synapse
%     Axis_z=[0 0 1];
%     norm_coeff=sqrt(SynR_coeff(1,3)^2+SynR_coeff(2,3)^2+SynR_coeff(3,3)^2); % normalize coeff
%     cos_a=Axis_z*SynR_coeff(:,3)/norm_coeff;
%     syn_alpha(rg)=acosd(abs(cos_a)); % arccos in degrees
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Find the center of synapse
%     radii=1;
%     [Ctr11,Sig11] = subclust(SynR_score,radii); % find center of a subcluster and the range of influence of the center in PCA coordinates
%     %%make sure Ctr0 has only 1 row
%     if length(Ctr11(:,1))==1
%         Ctr1=Ctr11;
%     else
%         Ctr1=mean(Ctr11);
%     end
%     
%     if length(Sig11(:,1)==1)
%         Sig1=Sig11;
%     else
%         Sig1=mean(Sig11);
%     end
%     Ctr0=mean(syn_rg,radii); %   find center of a subcluster and the range of influence of the center
%     %%make sure Ctr0 has only 1 row
%     %%plot cluster
% %     Sig0 = (radii .* (max(syn_rg) - min(syn_rg))) / sqrt(8.0);
%     scatter3(syn_rg(:,1), syn_rg(:,2), syn_rg(:,3),30,'filled','b'); % Plot scattered data in original space
%     
%     
%     % plot center
%     scatter3(Ctr0(1), Ctr0(2),Ctr0(3), 100, 'filled','r'); % Plot center of cluster
%     
%     % Pricipal axis
%     axis1_x1=-SynR_coeff(1,1)*500+Ctr0(1);
%     axis1_x2=SynR_coeff(1,1)*500+Ctr0(1);
%     axis1_y1=-SynR_coeff(2,1)*500+Ctr0(2);
%     axis1_y2=SynR_coeff(2,1)*500+Ctr0(2);
%     axis1_z1=-SynR_coeff(3,1)*500+Ctr0(3);
%     axis1_z2=SynR_coeff(3,1)*500+Ctr0(3);
%     plot3([axis1_x1 axis1_x2], [axis1_y1 axis1_y2], [axis1_z1 axis1_z2],'m','LineWidth',2)
%     % Second axis
%     axis2_x1=-SynR_coeff(1,2)*500+Ctr0(1);
%     axis2_x2=SynR_coeff(1,2)*500+Ctr0(1);
%     axis2_y1=-SynR_coeff(2,2)*500+Ctr0(2);
%     axis2_y2=SynR_coeff(2,2)*500+Ctr0(2);
%     axis2_z1=-SynR_coeff(3,2)*500+Ctr0(3);
%     axis2_z2=SynR_coeff(3,2)*500+Ctr0(3);
%     plot3([axis2_x1 axis2_x2], [axis2_y1 axis2_y2], [axis2_z1 axis2_z2],'c','LineWidth',2)
%     % Third axis, this is the main axis of a synapse, i.e. the synaptic
%     % orientation
%     axis3_x1=-SynR_coeff(1,3)*500+Ctr0(1);
%     axis3_x2=SynR_coeff(1,3)*500+Ctr0(1);
%     axis3_y1=-SynR_coeff(2,3)*500+Ctr0(2);
%     axis3_y2=SynR_coeff(2,3)*500+Ctr0(2);
%     axis3_z1=-SynR_coeff(3,3)*500+Ctr0(3);
%     axis3_z2=SynR_coeff(3,3)*500+Ctr0(3);
%     plot3([axis3_x1 axis3_x2], [axis3_y1 axis3_y2], [axis3_z1 axis3_z2],'k','LineWidth',2)
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % plot ellipsoid
%     [Ex,Ey,Ez] = ellipsoid(Ctr1(1),Ctr1(2),Ctr1(3),Sig1(1),Sig1(2),Sig1(3),20);
%     %surf(Ex,Ey,Ez)
%     %Reshape the grid for surface to make a matrix in X Y and Z for the
%     %ellipsoid
%     Ex_X=reshape(Ex,21*21,1);
%     Ey_Y=reshape(Ey,21*21,1);
%     Ez_Z=reshape(Ez,21*21,1);
%     EXYZ=[Ex_X,Ey_Y,Ez_Z];
%     % Rotate the ellipsoid back to in the original space
%     E_back=EXYZ*inv(SynR_coeff);
%     Ex_X_b=reshape(E_back(:,1),21,21);
%     Ey_Y_b=reshape(E_back(:,2),21,21);
%     Ez_Z_b=reshape(E_back(:,3),21,21);
%     surf(Ex_X_b,Ey_Y_b,Ez_Z_b)
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % create another ellipsoid with defined semi-axis such that within the
%     % ellisoid, the data is considered synaptic and outside of it, the data is
%     % considered extra-synaptic.
%     Syn_buff=1000; % The buffer layer outside synapse
%     [Syn_Ex, Syn_Ey, Syn_Ez]=ellipsoid(Ctr1(1),Ctr1(2),Ctr1(3),Sig1(1)+Syn_buff,Sig1(2)+Syn_buff,Sig1(3)+Syn_buff,20);
% %%%%%%%%%%%%%%%%%%%%%%%% dwell
% nearby_trace0=syn_nearby{rg};
% nTrace=length(nearby_trace0);
% dwell{rg}=[];
% for h=1:nTrace
%     dp=traces{nearby_trace0(h)}; % dimension of matrix of particle# K
% % nTrace=length(traces); %total particle number
% % dwell{rg}=[];
% % for h=1:nTrace
% %     dp=traces{h}; % dimension of matrix of particle# K
%     x1=dp(:,1);
%     y1=dp(:,2);
%     z1=dp(:,3);
%     frame1=dp(:,4);
%     inSyn=find((x1-Ctr0(1)).^2/(Sig1(1)+Syn_buff)^2+(y1-Ctr0(2)).^2/(Sig1(2)+Syn_buff)^2+(z1-Ctr0(3)).^2/(Sig1(3)+Syn_buff)^2<1);
%     head=1;
%     tail=1;
%     while head<length(inSyn)& tail<length(inSyn)
%         tail=tail+1;
%         if (inSyn(tail)-inSyn(tail-1))~=1
%             if tail-head>8
%                 dwell{rg}=vertcat(dwell{rg},[dp(1,5) frame1(tail-1)-frame1(head)]);
%                 head=tail;
%             end
%         end
%     end
%     if tail-head>8
%         dwell{rg}=vertcat(dwell{rg},[dp(1,5) frame1(tail-1)-frame1(head)]);
%     end
% end  
% %%%%%%%%%%%%%%%%%%%%%%%% dwell
%     
%     Syn_Ex_X=reshape(Syn_Ex,21*21,1);
%     Syn_Ey_Y=reshape(Syn_Ey,21*21,1);
%     Syn_Ez_Z=reshape(Syn_Ez,21*21,1);
%     Syn_EXYZ=[Syn_Ex_X,Syn_Ey_Y,Syn_Ez_Z];
%     % Rotate the ellipsoid back to in the original space
%     Syn_E_back=Syn_EXYZ*inv(SynR_coeff);
%     
%     Syn_Ex_b=reshape(Syn_E_back(:,1),21,21);
%     Syn_Ey_b=reshape(Syn_E_back(:,2),21,21);
%     Syn_Ez_b=reshape(Syn_E_back(:,3),21,21);
%     surf(Syn_Ex_b, Syn_Ey_b, Syn_Ez_b);
%     
%     alpha(.1)
%     
%     
%     %%In and out of synapse
%     clear k clear mTraces_rg traces2
%     % Get coordinate of traces near synapse_rg
%     if isempty(syn_nearby{rg})==0
%         for k=1:length(syn_trace_rg)
%             traces2{k}=traces{syn_trace_rg(k)};
%         end
%         mTraces_rg1=cell2mat(traces2');
%         mTraces_rg=mTraces_rg1(:,1:3);
%     else mTraces_rg1=zeros(5);mTraces_rg=mTraces_rg1(:,1:3);
%     end
%     
%     IN=inhull(mTraces_rg,Syn_E_back);
%     IN_ptc0=[IN, IN, IN].*mTraces_rg;
%     IN_ptc=IN_ptc0((IN_ptc0(:,1)~=0),:);
%     OUT_ptc0=mTraces_rg-IN_ptc0;
%     OUT_ptc=OUT_ptc0((OUT_ptc0(:,1)~=0),:);
%     scatter3(IN_ptc(:,1),IN_ptc(:,2),IN_ptc(:,3),30,'filled','r')
%     scatter3(OUT_ptc(:,1),OUT_ptc(:,2),OUT_ptc(:,3),30,'filled','g')
%     
%     
%     % Output
%     clear IN_trace0 IN_trace1
%     IN_trace0=[IN,IN,IN,IN,IN].*mTraces_rg1;
%     IN_trace1=IN_trace0((IN_trace0(:,1)~=0),:);
%     if isempty(IN_trace1)==0
%         IN_trace{rg}=IN_trace1;
%     else
%         IN_trace{rg}=zeros(1,5);
%     end
%     
%     OUT_trace0=mTraces_rg1-IN_trace0;
%     OUT_trace1=OUT_trace0((OUT_trace0(:,1)~=0),:);
%     if isempty(OUT_trace1)==0
%         OUT_trace{rg}=OUT_trace1;
%     else
%         OUT_trace{rg}=zeros(1,5);
%     end
%     
%     %%The line through synapse center and a point in the trace
%     IN_trace1_length=length(IN_trace1(:,1));
%     theta=zeros(1,IN_trace1_length); % make an array of zeros for theta
%     phi=zeros(1,IN_trace1_length);
%     for tr=1:IN_trace1_length
%         INtrace_point=[IN_trace1(tr,1)-Ctr0(1) IN_trace1(tr,2)-Ctr0(2) IN_trace1(tr,3)-Ctr0(3)]; % vector between the center of the synapse and a point in the trace
%         INtrace_X=INtrace_point*SynR_coeff(:,1);
%         INtrace_Y=INtrace_point*SynR_coeff(:,2);
%         INtrace_Z=INtrace_point*SynR_coeff(:,3);
%         
%         norm_coeff=sqrt(SynR_coeff(1,3)^2+SynR_coeff(2,3)^2+SynR_coeff(3,3)^2);
%         norm_INtrace=sqrt((IN_trace1(tr,1)-Ctr0(1))^2+(IN_trace1(tr,2)-Ctr0(2))^2+(IN_trace1(tr,3)-Ctr0(3))^2);
%         cos_theta=INtrace_Z/(norm_coeff*norm_INtrace);
%         sin_theta=sqrt(1-cos_theta^2);
%         norm_coeffX=sqrt(SynR_coeff(1,1)^2+SynR_coeff(2,1)^2+SynR_coeff(3,1)^2); % length of coeff at the 1st principal axis
%         cos_theta2=INtrace_X/(norm_coeffX*norm_INtrace); % angle between Intrace and the 1st pricipal axis
%         cos_phi=norm_INtrace*cos_theta2/(norm_INtrace*sin_theta);
%         theta(tr)=acosd(cos_theta); % theta: 0~pi
%         if INtrace_point*SynR_coeff(:,2)~=0
%             phi_0=acosd(cos_phi*sign(INtrace_Y)); % phi is only 0~ pi, actural phi should be 0~2pi, taking consideration of the sign on y axis
%         else phi_0=acosd(cos_phi);
%         end
%         
%         if INtrace_Y~=0
%             if INtrace_Y>0
%                 phi(tr)=phi_0;
%             else phi(tr)=phi_0+180;
%             end
%         else
%             phi(tr)=phi_0;
%         end
%         
%     end
%     %     figure
%      Omiga{rg}=[theta', phi'];
%     %     [pX,pY] = pol2cart(phi',theta');
%     %     compass(pX,pY)
%     
%     % % Test
%     %
%     % ptc=-500+1000*rand(100,3);
%     % ptc_ctr=[ptc(:,1)+Ctr0(1), ptc(:,2)+Ctr0(2), ptc(:,3)+Ctr0(3)]
%     % IN=inhull(ptc_ctr,Syn_E_back);
%     % IN_ptc0=[IN, IN, IN].*ptc_ctr;
%     % IN_ptc=IN_ptc0(find(IN_ptc0(:,1)),:);
%     % OUT_ptc0=ptc_ctr-IN_ptc0;
%     % OUT_ptc=OUT_ptc0(find(OUT_ptc0(:,1)),:);
%     %
%     % scatter3(IN_ptc(:,1),IN_ptc(:,2),IN_ptc(:,3),30,'filled','r')
%     % scatter3(OUT_ptc(:,1),OUT_ptc(:,2),OUT_ptc(:,3),30,'filled','g')
%     PHI{rg}=phi;
%     THETA{rg}=theta;
% end
% hold off
% saveas(gcf,strcat(userdirin,'f6'),'fig');
% for m=1:syn_n
%     figure(6+ceil(m/16.0))
%     subplot(4,4,m-ceil(m/16)*16+16)
%     polar(PHI{m}/360*2*pi,THETA{m}/360*2*pi,'.');
%     title(m);
% end
% 
%  save(strcat(userdirin,'dwell_time.mat'),'dwell');
%  save(strcat(userdirin,'syn_angle.mat'),'syn_alpha');
%   save(strcat(userdirin,'Omiga.mat'),'Omiga');
% clear Ctr0 Ctr00 Ctr1 Ctr11 EXYX E_back Ex Ex_X Ex_X_b Ey Ey_Y Ey_Y_b Ez Ez_Z Ez_Z_b
% clear Sig0 Sig00 Sig1 Sig11 SynR_coeff SynR_score Syn_EXYZ Syn_E_back Syn_Ex Syn_Ex_X Syn_Ex_b Syn_Ey Syn_Ey_Y Syn_Ey_b
% clear Syn_Ez Syn_Ez_Z Syn_Ez_b Syn_buff k latentmTraces_rg mTraces_rg1 radii rg syn_n syn_rg syn_rr
% clear axis1_x1 axis1_x2 axis1_y1 axis1_y2 axis1_z1 axis1_z2
% clear axis2_x1 axis2_x2 axis2_y1 axis2_y2 axis2_z1 axis2_z2
% clear axis3_x1 axis3_x2 axis3_y1 axis3_y2 axis3_z1 axis3_z2
% %% Step6
% %%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%trace_synapse_IN_OUT
% % This script output matrix of traces with give array of specific
% % particle ID for all traces near synapses and amoung these traces, the
% % traces (or part of traces that is) in the synapse and out of synapse.
% % This script need input from nearby_traces.m syn_COAR.m
% 
% %%Get trace ID
% % Must run nearby_trace.m before running this script.
% mTraceID=cell2mat(syn_nearby');
% TraceID=unique(mTraceID); % get the unique particle ID for the nearby traces and arrange in acsend order.
% 
% 
% 
% %%load file
% %make variable x y z frame number and particle ID
% x=result_tracking(:,1);
% y=result_tracking(:,2);
% z=result_tracking(:,3); % 0.79 is the correction factor for RI mismatch for water/     oil for an 1.40 NA objective
% frame=result_tracking(:,4);
% particle=result_tracking(:,5);
% d=length(particle); % total row number
% ID_length=length(TraceID); % length of particle number
% % Rearrange traces nearby according to ID and store it in Trace
% for k1=1:ID_length
%     k=TraceID(k1);
%     dp=find(particle==k); % dimension of matrix of particle# K
%     dplength=length(dp);
%     dpmax=max(dp); % for particleK, the first row number
%     dpmin=min(dp); % for particleK, the last row number
%     m(k)=k; % test k value
%     x1=x(dpmin:dpmax);
%     y1=y(dpmin:dpmax);
%     z1=z(dpmin:dpmax);
%     frame1=frame(dpmin:dpmax);
%     particle1=particle(dpmin:dpmax);
%     Trace{k}=[x1,y1,z1,frame1,particle1];
%     clear dp dplength dpmax dpmin x1 y1 z1 frame1
% end
% mTrace=cell2mat(Trace');
% mTraceIN0=cell2mat(IN_trace');
% mTraceIN1=mTraceIN0((mTraceIN0(:,5)~=0),:);
% mTraceIN=intersect(mTraceIN1, mTraceIN1,'rows','stable');
% mTraceOUT=setdiff(mTrace,mTraceIN,'rows','stable');
% save(strcat(userdirin,'mTraceIN.mat'),'mTraceIN');
% save(strcat(userdirin,'mTraceOUT.mat'),'mTraceOUT');
% clear x y z frame particle k d m x1 y1 z1 frame1 particle1 dp ID_length
% 
% %%%%%%%%%%%%%%%%
% %%%%%fMSD_m_IN
% % This script calculates the MSD of traces with give array of specific
% % particle ID.
% %%Get trace ID
% % Must run nearby_trace.m before running this script.
% 
% % mTraceIN0=cell2mat(IN_trace');
% % mTraceIN=mTraceIN0(find(mTraceIN0(:,5)),:);
% % TraceIN=unique(mTraceIN(:,5)); % get the unique particle ID for the nearby traces and arrange in acsend order.
% tic
% %%load file
% %make variable x y z frame number and particle ID
% x=mTraceIN(:,1);
% y=mTraceIN(:,2);
% z=mTraceIN(:,3); % 0.79 is the correction factor for RI mismatch for water/oil for an 1.40 NA objective
% frame=mTraceIN(:,4);
% particle=mTraceIN(:,5);
% d=length(particle); % total row number
% TraceIN=unique(mTraceIN(:,5));
% ID_length=length(TraceIN); % length of particle number
% tSteps=10;
% clear Dif
% %%MSD calculation
% figure
% for k1=1:ID_length
%     k=TraceIN(k1);
%     dp=find(particle==k); % dimension of matrix of particle# K
%     dplength=length(dp);
%     if dplength<2
%         Dif(k,:)=[0, 0];
%     else
%         dpmax=max(dp); % for particleK, the first row number
%         dpmin=min(dp); % for particleK, the last row number
%         m(k)=k; % test k value
%         x1=x(dpmin:dpmax);
%         y1=y(dpmin:dpmax);
%         z1=z(dpmin:dpmax);
%         frame1=frame(dpmin:dpmax);
%         particle1=particle(dpmin:dpmax);
%         Trace_IN{k}=[x1,y1,z1,frame1,particle1];
%         [MSD00,d2r0,counts]=fMSD_vect(x1,y1,z1,frame1,dpmax,dpmin,tSteps); % Call function MSD, MSD.m file must be in the same folder as this file
%         cutoff=10;  %cutoff for data points
%         ind=find(counts>cutoff-1); % find index of counts that is equate to and above cutoff
%         MSDCF=MSD00(ind); % find the MSD for those index from last line
%         %ind1=[0,ind]; % Add (0,0) as the first point of the curve
%         %MSDCF1=[0,MSDCF]; % Add (0,0) as the first point of the curve
%         % ind1=ind(1:4);
%         % MSDCF1=MSDCF(1:4);
%         indlength=length(ind);
%         if indlength>=4
%             ind1=ind(1:4);
%             MSDCF1=MSDCF(1:4);
%             Dif(k,:) = polyfit(ind1',MSDCF1',1);
% % ORIGINAL            f = fit(ind1',MSDCF1','poly1');
% % ORIGINAL            Dif(k,:)=coeffvalues(f); % First column is the diffusion coeff (slope). To get um^2/s multiply by 2*10^(-5), to get cm^2/s, multiply by 2*10^(-13)
%         else Dif(k,:)=[0,0]; ind1=ind; MSDCF1=MSDCF;
%         end
%     end
%     plot(ind1,MSDCF1)
%     hold all
%     clear dp dplength dpmax dpmin x1 y1 z1 frame1 particle1 MSD00 d2r0 counts  MSDCF f
% end
% hold off
% dt=0.05; %Default time step 0.05s=50ms
% Dif1=Dif(:,1);
% Dif_positive=Dif1(find(Dif1>0));
% Dif_time_IN=Dif_positive/(dt*2*3); %D for time step dt, for 3D
% Dif_track_IN=Dif_positive/(dt*2*3*10^6);
% figure
% hist(Dif_track_IN, 15);
% clear x y z frame particle k d m x1 y1 z1 frame1 dp Dif1 Dif_positive
% save(strcat(userdirin,'diffusion_coefficients_in.mat'),'Dif_track_IN');
% 
% toc
% fast_In_time=toc;
% %%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%fMSD_m_OUT
% % This script calculates the MSD of traces with give array of specific
% % particle ID.
% %%Get trace ID
% % Must run nearby_trace.m before running this script.
% % clear mTraceIN0 mTraceIN TraceIN
% % mTraceIN0=cell2mat(IN_trace');
% % mTraceIN=mTraceIN0(find(mTraceIN0(:,5)),:);
% % TraceIN=unique(mTraceIN(:,5)); % get the unique particle ID for the nearby traces and arrange in acsend order.
% % %%
% 
% %
% %This m_file seperates particles according to their ID in the 'result', and
% %calculates diffusion coeficient for each particle by fitting straight
% %lines to the first 4 points of MSD/dt of the particle.
% %of the tracking file and outputs the submatrices, e.g. sub_k is the
% %submatrix of particle# k.
% tSteps=10;
% %%load file
% %make variable x y z frame number and particle ID
% x=mTraceOUT(:,1);
% y=mTraceOUT(:,2);
% z=mTraceOUT(:,3); % 0.79 is the correction factor for RI mismatch for water/oil for an 1.40 NA objective
% frame=mTraceOUT(:,4);
% particle=mTraceOUT(:,5);
% d=length(particle); % total row number
% TraceOUT=unique(mTraceOUT(:,5));
% ID_length=length(TraceOUT); % length of particle number
% %%MSD calculation
% %%MSD calculation
% figure
% for k1=1:ID_length
%     k=TraceOUT(k1);
%     dp=find(particle==k); % dimension of matrix of particle# K
%     dplength=length(dp);
%     if dplength<2
%         Dif(k,:)=0;
%     else
%         dpmax=max(dp); % for particleK, the first row number
%         dpmin=min(dp); % for particleK, the last row number
%         m(k)=k; % test k value
%         x1=x(dpmin:dpmax);
%         y1=y(dpmin:dpmax);
%         z1=z(dpmin:dpmax);
%         frame1=frame(dpmin:dpmax);
%         particle1=particle(dpmin:dpmax);
%         Trace_OUT{k}=[x1,y1,z1,frame1,particle1];
%         % eval(['sub_' num2str(k) '=[x1,y1,z1,frame1]']);% creating submatrix for particle# k sub_k
%         [MSD00,d2r0,counts]=fMSD_vect(x1,y1,z1,frame1,dpmax,dpmin,tSteps); % Call function MSD, MSD.m file must be in the same folder as this file
%         cutoff=10;  %cutoff for data points
%         ind=find(counts>cutoff-1); % find index of counts that is equate to and above cutoff
%         MSDCF=MSD00(ind); % find the MSD for those index from last line
%         %ind1=[0,ind]; % Add (0,0) as the first point of the curve
%         %MSDCF1=[0,MSDCF]; % Add (0,0) as the first point of the curve
%         % ind1=ind(1:4);
%         % MSDCF1=MSDCF(1:4);
%         indlength=length(ind);
%         if indlength>=4
%             ind1=ind(1:4);
%             MSDCF1=MSDCF(1:4);
%             Dif(k,:) = polyfit(ind1',MSDCF1',1);
% % ORIGINAL            f = fit(ind1',MSDCF1','poly1');
% % ORIGINAL            Dif(k,:)=coeffvalues(f); % First column is the diffusion coeff (slope). To get um^2/s multiply by 2*10^(-5), to get cm^2/s, multiply by 2*10^(-13)        else Dif(k,:)=0;
%         end
%     end
%     plot(ind1,MSDCF1)
%     hold all
%     %eval(['MSD00_' num2str(k) '=[MSD00]']);
%     %eval(['d2r0_' num2str(k) '=[d2r0]']);
%     clear dp dplength dpmax dpmin x1 y1 z1 frame1
% end
% hold off
% dt=0.05; %Default time step 0.05s=50ms
% Dif1=Dif(:,1);
% Dif_positive=Dif1(find(Dif1>0));
% Dif_time_OUT=Dif_positive/(dt*2*3); %D for time step dt, for 3D
% Dif_track_OUT=Dif_positive/(dt*2*3*10^6);
% figure
% hist(Dif_track_OUT, 15);
% save(strcat(userdirin,'diffusion_coefficients_out.mat'),'Dif_track_OUT');

