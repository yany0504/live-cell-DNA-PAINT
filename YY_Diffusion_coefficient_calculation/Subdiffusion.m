%% this script is used to analyze diffusion of trajectories by fitting MSD curves to D*t^a
clc;%clear
close all;% Close all figures
if ~exist('fileName1','var')|| isempty(fileName1)
    [userfilein, userdirin]=uigetfile({
        '*.txt','Data file (*.xls)';...
        '*.*','All Files (*.*)'},'Select the Tracking file to process',...
        'E:\DATA_2\cell imaging\NMDA and AMPA receptors\NMDA receptor only\Live cell_tracking\20140918_NMDAR-705 nm+Homer1 mGeos\1\2');
    fileName1=fullfile(userdirin,userfilein);
else
    if ~exist(fileName1,'file')
        fprintf('File not found: %s\n',fileName1);
        return;
    end
end


result_tracking=textread(fileName1);
%%
dt=0.05; %Default time step 0.05s=50ms

tic
x=result_tracking(:,1);
y=result_tracking(:,2);
z=result_tracking(:,3); % 0.79 is the correction factor for RI mismatch for water/oil for an 1.40 NA objective
frame=result_tracking(:,4);
particle=result_tracking(:,5);
d=length(particle); %total row number
ptotal=particle(d);%total particle number
tSteps=100;
clear Dif
clear MSD_Ind
% figure
% hold all
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
    [MSD00,d2r0,counts,E_sem]=fMSD_vect_withError(x1,y1,z1,frame1,dpmax,dpmin,tSteps); % Call function MSD, MSD.m file must be in the same folder as this file
    cutoff=50;  %cutoff for data points
    ind=find(counts>cutoff-1); % find index of counts that is equate to and above cutoff
    MSDCF=MSD00(ind); % find the MSD for those index from last line
    ind_log=log10(ind*dt);
    MSDCF_log=log10(MSDCF*10^(-6));
    Err_sem=E_sem(ind);
%     figure
%     subplot(2,1,1)
%     plot(ind,MSDCF)
%     errorbar(MSDCF,Err_sem)
%     subplot(2,1,2)
%     plot(ind_log,MSDCF_log)
    indlength=length(ind);
    if indlength>=20
        ind1=ind_log(1:20);
        MSDCF1=MSDCF_log(1:20);
        subDif(j,:)= polyfit(ind1',MSDCF1',1);
        
% ORIGINAL        f = fit(ind1',MSDCF1','poly1');
% ORIGINAL       Dif(j,:)=coeffvalues(f); % First column is the diffusion coeff (slope). To get um^2/s multiply by 2*10^(-5), to get cm^2/s, multiply by 2*10^(-13)
    else Dif(j,:)=[0,0]; ind1=ind; MSDCF1=MSDCF; %%% Changed on 2013/05/13, Dif has to have 2 columns, and if indlenth<4, ind1 and MSDCF1 has to have values
    end
    particle2(j)=j; % Added by Anne, 1/14/2014
    plot(ind1,MSDCF1)
    
   
    inst_msd{j}= inst_MSD( x1,y1,z1,frame1);
end
% hold off
% saveas(gcf,strcat(userdirin,'MSD_sub'),'fig')

%%
% dt=0.05; %Default time step 0.05s=50ms
% Dif1=Dif(:,1);
% Dif_positive=Dif1(Dif1>0);
% Dif_time=Dif_positive/(dt*2*3); %D for time step dt, for 3D
% Dif_track=Dif_positive/(dt*2*3*10^6);
% Dif_track_ID=particle2(Dif1>0); % Added by Anne, 1/14/2014, for trace ID of corresponding D
% Dif_track_W_ID(:,1)=Dif_track;
% Dif_track_W_ID(:,2)=Dif_track_ID;
% figure
% hist(Dif_track, 30);
% saveas(gcf,strcat(userdirin,'Dif_track'),'fig')
subDif1(:,1)=nonzeros(subDif(:,1));
subDif1(:,2)=nonzeros(subDif(:,2));
save(strcat(userdirin,'subDif.mat'),'subDif1');
% save(strcat(userdirin,'diffusion_coefficients.mat'),'Dif_track_W_ID');
% save(strcat(userdirin,'Trace_range.mat'),'Trace_range')
% clear x y z frame particle j d m x1 y1 z1 frame1 dp
