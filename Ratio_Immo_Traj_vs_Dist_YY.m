% ratio of immobile trajectories 
clc;clear;close all;
FilePath = 'H:\Data_from07222021\DNA-PAINT(Live-cell)\SPT_tracking\Basal_tracking\SPT_results_0701_0708_0715\';
OutputFile='ND_DNA-PAINT';
criteria_slow=-2.5;
criteria_confined=-0.5;

dist_step=0:0.1:2;


%% load
if ~exist('FileName1','var')|| isempty(FileName1)
    [userfilein, userdirin]=uigetfile({
         '*.xlsx','Data file (*.xlsx)';...
        '*.*','All Files (*.*)'},'Select the diffusion vs trajectory file to process',...
        FilePath, 'MultiSelect','on');
    FileName1=fullfile(userdirin,userfilein);
else
    if ~exist(FileName1,'file')
        fprintf('File not found: %s\n',FileName1);
        return;
    end
end

if iscell(FileName1);
    for i=1:length(FileName1);
        data1{i}=xlsread(FileName1{i});
    end
else
    data1=xlsread(FileName1);
end

if iscell(data1);
    DiffTrajDist=data1{1};
    for i=2:length(FileName1);
        DiffTrajDist=cat(1,DiffTrajDist,data1{i});
    end
else
    DiffTrajDist=data1;
end

%%

if iscell(data1);
    ratio_slow=zeros(20,length(data1));
    num_data=length(data1);
    for j=1:num_data;
        
        for i=1:length(dist_step)-1;
            diff_range=data1{j}(data1{j}(:,5)>=dist_step(i)&data1{j}(:,5)<=dist_step(i+1),:);
            num_slow=sum(diff_range(:,2)<=criteria_slow);
            num_fast=sum(diff_range(:,2)>criteria_slow);
            num_total=num_slow+num_fast;
            frac_slow=num_slow/num_total;
            frac_fast=num_fast/num_total;
            ratio_slow(i,j)=frac_slow;
        end
        
    end
else
    ratio_slow=zeros(20,1);
    num_data=1;
        for i=1:length(dist_step)-1;
            diff_range=data1(data1(:,5)>=dist_step(i)&data1(:,5)<=dist_step(i+1),:);
            num_slow=sum(diff_range(:,2)<=criteria_slow);
            num_fast=sum(diff_range(:,2)>criteria_slow);
            num_total=num_slow+num_fast;
            frac_slow=num_slow/num_total;
            frac_fast=num_fast/num_total;
            ratio_slow(i,1)=frac_slow;
        end
    
end


frac_mean=mean(ratio_slow,2);
frac_std=std(ratio_slow,0,2);
frac_SEM=frac_std/length(data1);

x=0.1:0.1:2;
figure;
errorbar(x,frac_mean,frac_SEM,'k--','LineWidth',2);
hold on
scatter(x,frac_mean(:,1),'k','LineWidth',2);
set(gca,'FontSize',25,'fontweight','Bold');
xlabel('Distance from NN Homer1 (\mum)','FontSize',30);
ylabel('Fraction of GluA1 in nanodomains','FontSize',30);
set(gcf,'Position',[100 100 900 700])
saveas(gcf,[FilePath OutputFile],'fig');