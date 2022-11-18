%Calculate tau from in situ image by using dbscan
clc;clear;close all;
path='H:\Data_from07222021\DNA-PAINT(Live-cell)\Counting_qPAINT\cLTP_counting\cLTP_singlelabeling\210528_qPAINT_cLTP\Cell1\';
filename='GluA1_corrected.csv';
File=[path filename];
data=xlsread(File);
epsilon=10;
minpts=25;
exposure=0.1;
tic
idx=dbscan(data(:,3:4),epsilon,minpts);
toc


data_temp=data;
tic
for i=1:length(data);
    if idx(i)==-1;
        data_temp(i,:)=zeros(1,9);
    else
    end
end
toc
data_filtered=data_temp(data_temp(:,9)~=0,:);
idx_filtered=idx(idx~=-1);
figure;gscatter(data_filtered(:,3),data_filtered(:,4),idx_filtered);axis equal;
legend('off')
event_time=cell(max(idx_filtered),1);
for j=1:max(idx_filtered);
    event_time{j,1}=data_filtered(idx_filtered==j,:);
end
%%
DarkTime=zeros(length(event_time),2);
DarkTimeCI=zeros(length(event_time),2);
DarkTime(:,2)=1:length(DarkTime);
BrightTime=zeros(length(event_time),2);
BrightTimeCI=zeros(length(event_time),2);
BrightTime(:,2)=1:length(BrightTime);

for i=1:length(event_time);
    if size(event_time{i},1)<10 %isempty(events_synapse{i});
        fprintf('Not enough events for the synapse.\n')
    else
        [DarkTime(i,1) DarkTimeCI(i,:) BrightTime(i,1) BrightTimeCI(i,:)]=CalculateTauDark(event_time{i},exposure,5);
    end
  
end
DarkTime=DarkTime(DarkTime(:,1)~=0,:);
DarkTimeCI=DarkTimeCI(DarkTime(:,1)~=0,:);
BrightTime=BrightTime(BrightTime(:,1)~=0,:);
BrightTimeCI=BrightTimeCI(BrightTime(:,1)~=0,:);

%% test to get tau_Dark_1 by fitting GMM.
opts=statset('MaxIter',10000);% max number iterations
DarkTime1=DarkTime(DarkTime(:,1)<120&DarkTime(:,1)>2,:);
GMModel_DarkTime=fitgmdist(DarkTime1(:,1),2,'Options',opts);
X=0:1:150;
X=X';
y=pdf(GMModel_DarkTime,X);
% gmPDF_Dark=@(x) arrayfun(@(x0) pdf(GMModel_DarkTime,x0),x);

figure;h_dark=histogram(DarkTime1(:,1),'BinWidth',5,'Normalization','pdf');
hold on;
plot(X,y,'LineWidth',2);
set(gca,'FontSize',25,'LineWidth',2,'fontweight','Bold','YTick',[0,0.01,0.02])
xlabel('Dark time (s)','FontSize',30,'fontweight','Bold')
xlim([0 120])
ylabel('Frequency','FontSize',30,'fontweight','Bold')
% g2=gca;
% fplot(gmPDF_Dark,[g2.XLim],'LineWidth',3);
save([path filename(1:end-4) '_DarkTime=' num2str(max(GMModel_DarkTime.mu)) 'sec.mat']);
