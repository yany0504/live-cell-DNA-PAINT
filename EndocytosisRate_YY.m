clc;clear;close all;
FilePath = 'E:\Data\DNA-PAINT(Live-cell)\SPT_tracking\Basal_tracking\210624_live-cell_DNA-PAINT_func_of_time\Cell1\Local_endocytosis\';
OutputFile='test';

if ~exist('FileName1','var')|| isempty(FileName1)
    [userfilein, userdirin]=uigetfile({
         '*.csv','Data file (*.csv)';...
        '*.*','All Files (*.*)'},'Select dendrite file to process',...
        FilePath, 'MultiSelect','on');
    FileName1=fullfile(userdirin,userfilein);
else
    if ~exist(FileName1,'file')
        fprintf('File not found: %s\n',FileName1);
        return;
    else
    end
end

if iscell(FileName1);
    for i=1:length(FileName1);
        data1{i}=xlsread(FileName1{i});
    end
else
    data1=xlsread(FileName1);
end

if ~exist('FileName2','var')|| isempty(FileName2)
    [userfilein, userdirin]=uigetfile({
         '*.csv','Data file (*.csv)';...
        '*.*','All Files (*.*)'},'Select spine file to process',...
        FilePath, 'MultiSelect','on');
    FileName2=fullfile(userdirin,userfilein);
else
    if ~exist(FileName2,'file')
        fprintf('File not found: %s\n',FileName2);
        return;
    else
    end
end
if iscell(FileName2);
    for i=1:length(FileName2);
        data2{i}=xlsread(FileName2{i});
    end
else
    data2=xlsread(FileName2);
end

dendrite_frame=[];
spine_frame=[];

for i=1:length(data1);
    dendrite_frame=[dendrite_frame;data1{i}(:,2)];
end
for i=1:length(data2);
    spine_frame=[spine_frame;data2{i}(:,2)];
end

[muhat1,muci1]=mle(spine_frame,'dist','exp');

[muhat2,muci2]=mle(dendrite_frame,'dist','exp');
x=0:1:40000;
y_spine=1-exp(-x/muhat1);
y_dendrite=1-exp(-x/muhat2);
figure;
hs=histogram(spine_frame,'BinWidth',1000,'Normalization','cdf');
hold on;
plot(x,y_spine,'r')
xlabel('Frame');
ylabel('Cumulative Frequency');
title('Surface GluA1: spine');

figure;
hd=histogram(dendrite_frame,'BinWidth',1000,'Normalization','cdf');
hold on;
plot(x,y_dendrite,'r')
xlabel('Frame');
ylabel('Cumulative Frequency');
title('Surface GluA1: dendrite');

t_spine=0.05*muhat1;
conf_spine=0.05*muci1;
t_dendrite=0.05*muhat2;
conf_dendrite=0.05*muci2;