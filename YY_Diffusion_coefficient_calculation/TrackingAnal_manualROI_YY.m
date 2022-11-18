%% Step 0 Parameters. YOU NEED TO TAKE A LOOK AT THIS PART CAREFULLY, OTHERWISE YOU WILL GET ERRORS OR MISINTERPRETED RESULTS.
clc;clear;
% close all;

path = 'H:\Data_from07222021\DNA-PAINT(Live-cell)\SPT_tracking\cLTP_tracking\211009_liveDNA-PAINT_cLTP_noTRF\Cell2_vehicle\';
outputpath=[path '\Results'];
pathsplit=strsplit(path,'\');
outputindex=pathsplit(end-1);
outputname='before_mem3_good5_maxdisp500'; % this will go to all the output filex e.g. result_tracking_picasso.txt
isitpicasso=0;
actualdimension=3;
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
dt=0.05; %exposure time. Default time step 0.05s=50ms
tSteps=10; %length of MSD
diffusion_thres = 4; % how many time points used in calculating diffusion coefficients by linear fitting MSD
cutoff=param.good; % in MSD calculation minimum length of trajectory line 249
%for cluster analysis
% k - number of objects in a neighborhood of an object
% (minimal number of objects considered as a cluster)
% Eps - neighborhood radius, if not known avoid this parameter or put []
filter_photon=0; % if you filter non-specific binding using photon number 1, otherwise 0.

%% load file
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

data_QD=xlsread(fileName2);
data_QD(isnan(data_QD))=Inf;
if actualdimension==3;
    data_QD=data_QD(abs(data_QD(:,5))<threshold,:);
    data_QD(:,5)=data_QD(:,5)*factor;
else actualdimension==2;
end
data_QD=sortrows(data_QD,2);



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

% param.mem=80;
% param.dim=3; %for 3D use 3
% param.good=10; %original 50
% param.quiet=0;
% maxdisp=500;

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
%     Dif_time=Dif_positive/(dt*2*3); %D for time step dt, for 3D
    Dif_track=Dif_positive/(dt*2*3*10^6);
else actualdimension==2;
%         Dif_time=Dif_positive/(dt*2*2); %D for time step dt, for 2D
    Dif_track=Dif_positive/(dt*2*2*10^6);
end
Dif_track_ID=particle2(Dif1>0); % Added by Anne, 1/14/2014, for trace ID of corresponding D
Dif_track_W_ID(:,1)=Dif_track;
Dif_track_W_ID(:,2)=Dif_track_ID;

figure;
hDiff=histogram(log10(Dif_track),'NumBins',30,'Normalization','pdf');
xlim([-5 0])
xlabel('Diffusion coefficient (\mum^2/s, log)');
ylabel('Frequency');
%hist(log10(Dif_track), 30);
saveas(gcf,strcat(outputpath,['Dif_track_' outputname]),'fig')
save(strcat(outputpath,['Inst_diffusion_coefficients_' outputname '.mat']),'inst_msd');
save(strcat(outputpath,['diffusion_coefficients_' outputname '.mat']),'Dif_track_W_ID');
%YY commented save(strcat(userdirin,'Trace_range.mat'),'Trace_range')
clear x y z frame particle j d m x1 y1 z1 frame1 dp