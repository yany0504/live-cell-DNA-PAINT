% Calculate counting error vs. frame number
clc;clear;close all;
FilePath = 'E:\Data\DNA-PAINT(Live-cell)\cLTP_counting\201130_livecellDNA-PAINT_STP\1\Results\';
OutputFile='Count_before';
if ~exist('FileName1','var')|| isempty(FileName1)
    [userfilein, userdirin]=uigetfile({
         '*.csv','Data file (*.csv)';...
        '*.*','All Files (*.*)'},'Select the diffusion vs trajectory file to process',...
        FilePath, 'MultiSelect','on');
    FileName1=fullfile(userdirin,userfilein);
else
    if ~exist(FileName1,'file')
        fprintf('File not found: %s\n',FileName1);
        return;
    end
end
Data=xlsread([userdirin userfilein]);
%binsize=[100:100:2000]';
binsize=[100 200 300 400 500 1000 2000 2500 4000]';
events=cell(length(binsize),1);
for n=1:length(binsize);
    [events{n} edge]=histcounts(Data(:,2),'BinWidth',binsize(n,1));
end

mean_binsize=zeros(length(events),1);
std_binsize=zeros(length(events),1);
for i=1:length(binsize);
    mean_binsize(i)=mean(events{i});
    std_binsize(i)=std(events{i});
end
percent_error=100*(std_binsize./mean_binsize);

figure;
scatter(binsize,percent_error);
ax=gca;
ax.FontSize=20;
xlabel('Bin size (frame)','FontSize', 30);
ylabel('% Error','FontSize', 30);
title(OutputFile);
set(gcf,'Position',[100, 100, 1000, 800]);
