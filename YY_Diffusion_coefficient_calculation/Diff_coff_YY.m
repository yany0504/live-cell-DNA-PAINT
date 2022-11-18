%% Getting Diffusion coefficients
%% Step1
clc;
clear
close all;% Close all figures

%% Reading QD file
path = 'G:\Data\BrainSlices\GluR2_antibody\181005_overexpressiontest_neurons\4\6';
if ~exist('fileName2','var')|| isempty(fileName2)
    [userfilein, userdirin]=uigetfile({
        '*.txt','Data file (*.txt)';...
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

if ~exist('threshold','var')||isempty(threshold)||~exist('factor','var')||isempty(factor)
    prompt = {'Enter threshold:','Enter factor of z:'};
    dlg_title = 'Input for preparing data';
    num_lines = 1;
    defaultanswer={'600','0.79'};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultanswer);
    threshold=str2double(answer(1));
    factor=str2double(answer(2));
end

data_QD=dlmread(fileName2);
data_QD(isnan(data_QD))=Inf;
data_QD=data_QD(abs(data_QD(:,7))<threshold,:);
data_QD(:,7)=data_QD(:,7)*factor;
data_QD=sortrows(data_QD,15);

%% Setting params
x=data_QD(:,5);
y=data_QD(:,6);
z=data_QD(:,7);
t=data_QD(:,15);
positionlist=[x y z t];
param.mem=10; %original was 10
param.dim=3; %dimensionality
param.good=50; %original was 50
param.quiet=0; %not affect result
maxdisp=500;
%% tracking
result_tracking = track( positionlist, maxdisp, param );


fid=fopen([userdirin 'result_tracking.txt'],'w');
fprintf(fid,'%f %f %f %d %d\n',result_tracking');    
fclose(fid);
%ORIGINAL save('result_tracking','-ascii');
clear x y z t positionlist param maxdisp

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
%saveas(gcf,strcat(userdirin,'MSD'),'fig')

dt=0.05; %Default time step 0.05s=50ms
Dif1=Dif(:,1);
Dif_positive=Dif1(Dif1>0);
Dif_time=Dif_positive/(dt*2*3); %D for time step dt, for 3D
Dif_track=Dif_positive/(dt*2*3*10^6);
Dif_track_ID=particle2(Dif1>0); % Added by Anne, 1/14/2014, for trace ID of corresponding D
Dif_track_W_ID(:,1)=Dif_track;
Dif_track_W_ID(:,2)=Dif_track_ID;
figure
hist(Dif_track, 30);
saveas(gcf,strcat(userdirin,'Dif_track'),'fig')
save(strcat(userdirin,'Inst_diffusion_coefficients.mat'),'inst_msd');
save(strcat(userdirin,'diffusion_coefficients.mat'),'Dif_track_W_ID');
save(strcat(userdirin,'Trace_range.mat'),'Trace_range')
clear x y z frame particle j d m x1 y1 z1 frame1 dp
