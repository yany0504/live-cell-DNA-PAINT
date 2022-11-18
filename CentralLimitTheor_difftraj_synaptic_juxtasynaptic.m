% Central limit Theorem: for large n, sample mean distribution is normal.

%% load files
clc;clear;close all;
FilePath = 'E:\Data\AMPARtracking\210226_0211and0224_SPT-QD_results\';
OutputFile='scFv-QD_CLT';
intra=0.5;
juxta=2;

if ~exist('FileName1','var')|| isempty(FileName1)
    [userfilein, userdirin]=uigetfile({
         '*.xlsx','Data file (*.xlsx)';...
        '*.*','All Files (*.*)'},'Select the diffusion, trajectory, distance file to process',...
        FilePath, 'MultiSelect','on');
    FileName1=fullfile(userdirin,userfilein);
else
    if ~exist(FileName1,'file')
        fprintf('File not found: %s\n',FileName1);
        return;
    end
end

for i=1:length(FileName1);
    data1{i,1}=xlsread(FileName1{i});
end
%% Sort synaptic, juxtasynaptic
synaptic=cell(length(data1),1);
juxtasynaptic=cell(length(data1),1);
extrasynaptic=cell(length(data1),1);
for i=1:length(data1);
    synaptic{i}=data1{i}(data1{i}(:,5)<intra,:);
    juxtasynaptic{i}=data1{i}(data1{i}(:,5)>intra&data1{i}(:,5)<juxta,:);
    extrasynaptic{i}=data1{i}(data1{i}(:,5)>juxta,:);
end

synaptic_means=zeros(length(data1),4);
juxtasynaptic_means=zeros(length(data1),4);

for i=1:length(data1);
    synaptic_means(i,1)=mean(synaptic{i}(:,2));
    synaptic_means(i,2)=std(synaptic{i}(:,2));
    synaptic_means(i,3)=mean(synaptic{i}(:,4));
    synaptic_means(i,4)=std(synaptic{i}(:,4));
    juxtasynaptic_means(i,1)=mean(juxtasynaptic{i}(:,2));
    juxtasynaptic_means(i,2)=std(juxtasynaptic{i}(:,2));
    juxtasynaptic_means(i,3)=mean(juxtasynaptic{i}(:,4));
    juxtasynaptic_means(i,4)=std(juxtasynaptic{i}(:,4));
end
%% plot
close all;
[p_diff,t_diff,stats_diff]=anova1([synaptic_means(:,1), juxtasynaptic_means(:,1)]);
xticklabels({'Synaptic', 'Juxtasynaptic'});
ylabel('Diffusion coefficient (\mum^2/s, log)');
set(gca,'FontSize',15,'fontweight','bold');
set(gcf, 'Position', [100 100 650 500]);
DataX = interp1( [0 1], xlim(), 0.1 );
DataY = interp1( [0 1], ylim(), 0.8 );
text(DataX,DataY,['p-value = ' num2str(round(p_diff,4))],'FontSize',15);
saveas(gcf,[FilePath OutputFile '_ANOVA_diffusion'],'fig');

[p_traj,t_traj,stats_traj]=anova1([synaptic_means(:,3), juxtasynaptic_means(:,3)]);
xticklabels({'Synaptic', 'Juxtasynaptic'});
ylabel('Trajectory range (\mum, log)');
set(gca,'FontSize',15,'fontweight','bold');
set(gcf, 'Position', [100 100 650 500]);
DataX = interp1( [0 1], xlim(), 0.1 );
DataY = interp1( [0 1], ylim(), 0.8 );
text(DataX,DataY,['p-value = ' num2str(round(p_traj,4))],'FontSize',15);
saveas(gcf,[FilePath OutputFile '_ANOVA_trajectory'],'fig');

%%
diff_syn=[mean(synaptic_means(:,1)),std(synaptic_means(:,1))];
diff_jux=[mean(juxtasynaptic_means(:,1)),std(juxtasynaptic_means(:,1))];
traj_syn=[mean(synaptic_means(:,3)),std(synaptic_means(:,3))];
traj_jux=[mean(juxtasynaptic_means(:,3)),std(juxtasynaptic_means(:,3))];