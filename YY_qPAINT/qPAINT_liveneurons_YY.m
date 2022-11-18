% qPAINT for live neurons with result_tracking and Diff_Traj_Dist files
%% load files
clear;clc;close all;

FilePath='C:\Users\yyoun4\Desktop\20210416_YL_scFv_w_5xR1_500pM_R1-9nt-Cy3B_100ms';

if ~exist('FileName_Diff_traj_dist','var')|| isempty(FileName_Diff_traj_dist)
    [userfilein1, userdirin1]=uigetfile({
         '*.xlsx','Data file (*.xlsx)';...
        '*.*','All Files (*.*)'},'Select Diff_Traj_Dist file to process',...
        FilePath);
    FileName_Diff_traj_dist=fullfile(userdirin1,userfilein1);
else
    if ~exist(FileName_Diff_traj_dist,'file')
        fprintf('File not found: %s\n',FileName_Diff_traj_dist);
        return;
    end
end

if ~exist('FileName_resulttracking','var')|| isempty(FileName_resulttracking)
    [userfilein2, userdirin2]=uigetfile({
         '*.txt','Data file (*.txt)';...
        '*.*','All Files (*.*)'},'Select result_tracking file to process',...
        FilePath);
    FileName_resulttracking=fullfile(userdirin2,userfilein2);
else
    if ~exist(FileName_resulttracking,'file')
        fprintf('File not found: %s\n',FileName_resulttracking);
        return;
    end
end

difftrajdist=xlsread([userdirin1 userfilein1]);
difftrajdist=difftrajdist(difftrajdist(:,5)<0.5,:);
tracking=textread([userdirin2 userfilein2]);

%% target synapses and trajectories
histogram(difftrajdist(:,8));
mobile_criteria=-1.8;
lengthlimit=100;
num_syn=23;
target_Traj=difftrajdist(find(difftrajdist(:,8)==num_syn),:);
target_Traj=target_Traj(target_Traj(:,7)<lengthlimit,:);
target_mobile=target_Traj(target_Traj(:,2)>mobile_criteria,:);
target_immobile=target_Traj(target_Traj(:,2)<=mobile_criteria,:);
traj_mobile=cell(size(target_mobile,1),1);
for i=1:length(traj_mobile);
    traj_mobile{i}=tracking(find(tracking(:,5)==target_mobile(i,6)),:);
end

traj_immobile=cell(size(target_immobile,1),1);
for i=1:length(traj_immobile);
    traj_immobile{i}=tracking(find(tracking(:,5)==target_immobile(i,6)),:);
end

%% calculate tau bright and tay dark for each fraction.
exposure=0.05;
frame=2000;
time=linspace(1,frame,frame)'*exposure;
intensity_mobile=zeros(frame,1);
for i=1:length(traj_mobile);
    for j=traj_mobile{i}(1,4):traj_mobile{i}(end,4);
        intensity_mobile(j)=intensity_mobile(j)+1;
 
    end
end

intensity_immobile=zeros(frame,1);
for i=1:length(traj_immobile);
    for j=traj_immobile{i}(1,4):traj_immobile{i}(end,4);
        intensity_immobile(j)=intensity_immobile(j)+1;
 
    end
end
figure;plot(time,intensity_mobile);
figure;plot(time,intensity_immobile);
[DarkTime_mobile,Darkstd_mobile, BrightTime_mobile,Brightstd_mobile]=calc_tau(intensity_mobile,exposure);
[DarkTime_immobile,Darkstd_immobile, BrightTime_immobile,Brightstd_immobile]=calc_tau(intensity_immobile,exposure);

