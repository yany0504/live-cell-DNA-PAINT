% visualize nanodomains based on GM model derived from cluster analysis
% "Distancce_2Dcontour_MLEverison_xxxxxx.m".
% input: result tracking file for a cell, synapse file, diffusion
% coefficient file, MLE_xxx_GMModel_analysis.mat file.

%% Step 1: location label....
clc;clear;close all;
FilePath = 'H:\Data_from07222021\DNA-PAINT(Live-cell)\SPT_tracking\Basal_tracking\210507_newscFv-5xR1_test_HP_DIV16\1\';
OutputFile='testingND\testresults';

%% load files.
if ~exist('FileName_GMModel','var')|| isempty(FileName_GMModel)
    [userfilein_GMModel, userdirin_GMModel]=uigetfile({
         '*.mat','Data file (*.mat)';...
        '*.*','All Files (*.*)'},'Select Gaussian Mixture Model results to process',...
        FilePath, 'MultiSelect','on');
    FileName_GMModel=fullfile(userdirin_GMModel,userfilein_GMModel);
else
    if ~exist(FileName_GMModel,'file')
        fprintf('File not found: %s\n',FileName_GMModel);
        return;
    end
end
load([FileName_GMModel],'GMModel_syn');
if ~exist('FileName_Track','var')|| isempty(FileName_Track)
    [userfilein_Track, userdirin_Track]=uigetfile({
         '*.txt','Data file (*.txt)';...
        '*.*','All Files (*.*)'},'Select result_tracking file to process',...
        FilePath, 'MultiSelect','on');
    FileName_Track=fullfile(userdirin_Track,userfilein_Track);
else
    if ~exist(FileName_Track,'file')
        fprintf('File not found: %s\n',FileName_Track);
        return;
    end
end
result_tracking=textread(FileName_Track);
if ~exist('FileName_Synapse','var')|| isempty(FileName_Synapse)
    [userfilein_Synapse, userdirin_Synapse]=uigetfile({
         '*.mat','Data file (*.mat)';...
        '*.*','All Files (*.*)'},'Select Synapse1 file to process',...
        FilePath, 'MultiSelect','on');
    FileName_Synapse=fullfile(userdirin_Synapse,userfilein_Synapse);
else
    if ~exist(FileName_Synapse,'file')
        fprintf('File not found: %s\n',FileName_Synapse);
        return;
    end
end
load(FileName_Synapse);
if ~exist('FileName_traceR','var')|| isempty(FileName_traceR)
    [userfilein_traceR, userdirin_traceR]=uigetfile({
         '*.mat','Data file (*.mat)';...
        '*.*','All Files (*.*)'},'Select trace_R file to process',...
        FilePath, 'MultiSelect','on');
    FileName_traceR=fullfile(userdirin_traceR,userfilein_traceR);
else
    if ~exist(FileName_traceR,'file')
        fprintf('File not found: %s\n',FileName_traceR);
        return;
    end
end
load(FileName_traceR);

if ~exist('FileName_DiffTrajDist','var')|| isempty(FileName_DiffTrajDist)
    [userfilein_DiffTrajDist, userdirin_DiffTrajDist]=uigetfile({
         '*.xlsx','Data file (*.xlsx)';...
        '*.*','All Files (*.*)'},'Select DiffTrajDist to process',...
        FilePath, 'MultiSelect','on');
    FileName_DiffTrajDist=fullfile(userdirin_DiffTrajDist,userfilein_DiffTrajDist);
else
    if ~exist(FileName_DiffTrajDist,'file')
        fprintf('File not found: %s\n',FileName_DiffTrajDist);
        return;
    end
end
DiffTrajDist=xlsread(FileName_DiffTrajDist);


%% total: clustering using GMModel
grp=cluster(GMModel_syn,[DiffTrajDist(:,2),DiffTrajDist(:,4)]);
figure;
gscatter(DiffTrajDist(:,2),DiffTrajDist(:,4),grp,'bgm');
ylabel('Trajectory range (log, \mum)');
xlabel('Diffusion coefficient (log, /mum^2/s)');
title('Total GMM cell1 06/10/2021');
xlim([-6 0]);ylim([-2 0.5]);
%% total: take group 2.
DTD_immo=DiffTrajDist(grp==2,:);
figure;
scatter(DTD_immo(:,2),DTD_immo(:,4),'g','filled');
ylabel('Trajectory range (log, \mum)');
xlabel('Diffusion coefficient (log, /mum^2/s)');
title('Immobile cell1 06/10/2021');
xlim([-6 0]);ylim([-2 0.5]);
%% filter synaptic = distance =< 0.5 um
DiffTrajDist_syn=DiffTrajDist(DiffTrajDist(:,5)<=0.5,:);


%% synaptic: clustering using GMModel
grp_syn=cluster(GMModel_syn,[DiffTrajDist_syn(:,2),DiffTrajDist_syn(:,4)]);
figure;
gscatter(DiffTrajDist_syn(:,2),DiffTrajDist_syn(:,4),grp_syn,'bgm');
ylabel('Trajectory range (log, \mum)');
xlabel('Diffusion coefficient (log, /mum^2/s)');
title('Synaptic GMM cell1 06/10/2021');
xlim([-6 0]);ylim([-2 0.5]);

%% synaptic: take group 2.
DTD_ND=DiffTrajDist_syn(grp_syn==2,:);
figure;
scatter(DTD_ND(:,2),DTD_ND(:,4),'g','filled');
ylabel('Trajectory range (log, \mum)');
xlabel('Diffusion coefficient (log, /mum^2/s)');
title('Synaptic & immobile cell1 06/10/2021');
xlim([-6 0]);ylim([-2 0.5]);
%% # nanodomains/synapse---
ID_synapse=unique(DTD_ND(:,8));
synapses=[1:length(Synapses1)]';
count=histc(DTD_ND(:,8),synapses);
figure;
h=histogram(count); % number of nanodomains / synapse histogram


%% Step 2. generate colormap for diffusion
Diffmin=-5;
Diffmax=-0.5;
Diffstep=10;
figure;
ax = axes;
cmap=colormap(jet(Diffstep));
c=colorbar('Ticks',[0,1],'TickLabels',{['10^' '{' num2str(Diffmin) '}'], ['10^' '{' num2str(Diffmax) '}']},'FontSize',20,'Location','west');
c.Label.String='Diffusion coefficient (\mum^2/s)';
ax.Visible = 'off';
Diff_range=linspace(Diffmin+abs(Diffmin+Diffmax)/(2*Diffstep),Diffmax-abs(Diffmin+Diffmax)/(2*Diffstep),Diffstep);
