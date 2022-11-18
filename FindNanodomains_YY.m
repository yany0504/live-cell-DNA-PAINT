% code to find nanodomains and its ROI (1x1 um2 window around Homer center)
% Input: Synapse1 file -> list of ROI
% Input2: SMLM of GluA1 and Homer -> visualize them.
%% location of file
clc;clear;%close all;
FilePath = 'E:\Data\AMPARtracking\210121_SPT_GluA1scFv\1to1QD\4\';
OutputFile='4QD';
box_size=500;
timefilter=0;
time_i=1;
time_l=500;

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
Synapses1_temp=Synapses1;
% localize synapses for generating ROIs.
for i=1:length(Synapses1_temp);
    center_x(i,1)=mean(Synapses1_temp{i}(:,1));
    center_y(i,1)=mean(Synapses1_temp{i}(:,2));
    center_z(i,1)=mean(Synapses1_temp{i}(:,3));
end
center_ROI=[center_x, center_y,center_z];

% load of ThunderSTORM image.
if ~exist('fileName1','var')|| isempty(fileName1)
    [userfilein1, userdirin1]=uigetfile({
         '*.csv','Data file (*.csv)';...
        '*.*','All Files (*.*)'},'Select the Homer file to process',...
        FilePath);
    fileName1=fullfile(userdirin1,userfilein1);
else
    if ~exist(fileName1,'file')
        fprintf('File not found: %s\n',fileName1);
        return;
    end
end
Homer=xlsread(fileName1);


if ~exist('fileName2','var')|| isempty(fileName2)
    [userfilein2, userdirin2]=uigetfile({
         '*.csv','Data file (*.csv)';...
        '*.*','All Files (*.*)'},'Select the DNA-PAINT file to process',...
        FilePath);
    fileName2=fullfile(userdirin2,userfilein2);
else
    if ~exist(fileName2,'file')
        fprintf('File not found: %s\n',fileName2);
        return;
    else [userdirin2,~,~]=fileparts(fileName2);
        userdirin2=strcat(userdirin2,'\');
    end
end
AMPAR1=xlsread(fileName2);
if timefilter==1;
    AMPAR=AMPAR1(AMPAR1(:,2)>=time_i&AMPAR1(:,2)<=time_l,:);
else
    AMPAR=AMPAR1;
end
%% plot synapses and draw ROI
figure
for n=1:length(Synapses1_temp)
        scatter3(Synapses1_temp{n}(:,1),Synapses1_temp{n}(:,2),Synapses1_temp{n}(:,3),20,'filled');
        a=sprintf('%d',n);
        text(center_ROI(n,1)+100,center_ROI(n,2),center_ROI(n,3),a);
        hold all
end
hold off
axis equal
view(-360,-90);
roi=images.roi.Polygon;
draw(roi)

%% filter synapses in ROI.
in_roi=inpolygon(center_ROI(:,1),center_ROI(:,2),roi.Position(:,1),roi.Position(:,2));

for i=1:length(Synapses1_temp);
    if in_roi(i)==1
    else
        Synapses1_temp{i}=[];
        center_ROI(i,:)=[0,0,0];
    end
end
Synapses1_temp=Synapses1_temp(~cellfun('isempty',Synapses1_temp));

%% localize synapses for generating ROIs.
for i=1:length(Synapses1_temp);
    center_x_new(i,1)=mean(Synapses1_temp{i}(:,1));
    center_y_new(i,1)=mean(Synapses1_temp{i}(:,2));
    center_z_new(i,1)=mean(Synapses1_temp{i}(:,3));
end
center_ROI_new=[center_x_new, center_y_new,center_z_new];

figure;
for n=1:length(Synapses1_temp)
        scatter3(Synapses1_temp{n}(:,1),Synapses1_temp{n}(:,2),Synapses1_temp{n}(:,3),20,'filled');
        a=sprintf('%d',n);
        text(center_ROI_new(n,1)+100,center_ROI_new(n,2),center_ROI_new(n,3),a);
        hold all
end
hold off
axis equal
view(-360,-90);

%% collect localizations of Homer and AMPARs for each ROI.
Homer_ROI=cell(length(center_ROI_new),1);
AMPAR_ROI=cell(length(center_ROI_new),1);

for i=1:length(center_ROI_new);
    AMPAR_ROI{i}=AMPAR(AMPAR(:,3)>=center_ROI_new(i,1)-box_size&AMPAR(:,3)<=center_ROI_new(i,1)+box_size&AMPAR(:,4)>=center_ROI_new(i,2)-box_size&AMPAR(:,4)<=center_ROI_new(i,2)+box_size,:);
    Homer_ROI{i}=Homer(Homer(:,3)>=center_ROI_new(i,1)-box_size&Homer(:,3)<=center_ROI_new(i,1)+box_size&Homer(:,4)>=center_ROI_new(i,2)-box_size&Homer(:,4)<=center_ROI_new(i,2)+box_size,:);
end


%% statistics for entire synapses 2D
% epsilon=75; % radius for cluster
% minpts=30; % minimum number of points for core points
% opts=statset('MaxIter',10000);% max number iterations
% 
% number_ND=zeros(length(AMPAR_ROI),1); 
% long_axis=[];
% short_axis=[];
% for i=1:length(AMPAR_ROI);
%     try
%         [idx,corepts]=dbscan(AMPAR_ROI{i}(:,3:4),epsilon,minpts);
%         data_temp=AMPAR_ROI{i};
%         for j=1:size(data_temp,1);
%             if idx(j)==-1;
%                 data_temp(j,:)=zeros(1,12);
%             else
%             end
%         end
%         data_filtered=data_temp(data_temp(:,9)~=0,:);
%         idx_filtered=idx(idx~=-1);
%         num_peaks=length(unique(idx_filtered));
%         gmmodel=fitgmdist([data_filtered(:,3),data_filtered(:,4)],num_peaks,'Options',opts);
%         size_ND=zeros(num_peaks,2);
%         for k=1:num_peaks;
%             size_ND(k,1)=2.3*sqrt(gmmodel.Sigma(1,1,k));
%             size_ND(k,2)=2.3*sqrt(gmmodel.Sigma(2,2,k));
%         end
%         long_axis1=max(size_ND')';
%         short_axis1=min(size_ND')';
%         long_axis=[long_axis;long_axis1];
%         short_axis=[short_axis;short_axis1];
%         % quantification of number
%         number_ND(i)=length(unique(idx_filtered));
%     catch
%         warning('No nanodomains found. Assign a value of 0.');
%     end
% end
% figure;
% h=histogram(number_ND,'BinWidth',1,'Normalization','probability');
% xlim([0 8])
% set(gca,'FontSize',20,'FontWeight','Bold')
% xlabel('# of Nanodomains','FontSize',30)
% ylabel('Frequency','FontSize',30)
% set(gcf,'position',[100 100 800 600]);
% 
% figure;
% h1=histogram(long_axis,'BinWidth',50,'Normalization','probability');
% xlim([0 600])
% set(gca,'FontSize',20,'FontWeight','Bold')
% xlabel('Length of principal axis (nm)','FontSize',30)
% %ylabel('Frequency','FontSize',30)
% set(gcf,'position',[100 100 800 600])
% 
% figure;
% h2=histogram(short_axis,'BinWidth',50,'Normalization','probability');
% xlim([0 600])
% set(gca,'FontSize',20,'FontWeight','Bold')
% xlabel('Length of auxiliary axis (nm)','FontSize',30)
% %ylabel('Frequency','FontSize',30)
% set(gcf,'position',[100 100 800 600])
% 
% save(strcat(FilePath, 'Nanodomain_', OutputFile, '.mat'));

%% statistics for entire synapses 3D
epsilon=100; % radius for cluster
minpts=25; % minimum number of points for core points
opts=statset('MaxIter',10000);% max number iterations

number_ND=zeros(length(AMPAR_ROI),1); 
long_axis=[];
median_axis=[];
short_axis=[];
for i=1:length(AMPAR_ROI);
    try
        if isempty(AMPAR_ROI{i})
        else
            [idx,corepts]=dbscan(AMPAR_ROI{i}(:,3:5),epsilon,minpts);
            data_temp=AMPAR_ROI{i};
            for j=1:size(data_temp,1);
                if idx(j)==-1;
                    data_temp(j,:)=zeros(1,12);
                else
                end
            end
            data_filtered=data_temp(data_temp(:,9)~=0,:);
            idx_filtered=idx(idx~=-1);
            num_peaks=length(unique(idx_filtered));
            if num_peaks == 0;
            else
                gmmodel=fitgmdist([data_filtered(:,3),data_filtered(:,4),data_filtered(:,5)],num_peaks,'Options',opts);
                size_ND=zeros(num_peaks,3);
                for k=1:num_peaks;
                    size_ND(k,1)=2.3*sqrt(gmmodel.Sigma(1,1,k));
                    size_ND(k,2)=2.3*sqrt(gmmodel.Sigma(2,2,k));
                    size_ND(k,3)=2.3*sqrt(gmmodel.Sigma(3,3,k));
                end
                long_axis1=max(size_ND')';
                median_axis1=median(size_ND')';
                short_axis1=min(size_ND')';
                long_axis=[long_axis;long_axis1];
                median_axis=[median_axis;median_axis1];
                short_axis=[short_axis;short_axis1];
                % quantification of number
                number_ND(i)=length(unique(idx_filtered));
            end
            
        end
    catch
        warning('No nanodomains found. Assign a value of 0.');
    end
end
figure;
h=histogram(number_ND,'BinWidth',1,'Normalization','probability');
xlim([0 8])
set(gca,'FontSize',20,'FontWeight','Bold')
xlabel('# of Nanodomains','FontSize',30)
ylabel('Frequency','FontSize',30)
set(gcf,'position',[100 100 800 600]);

figure;
h1=histogram(long_axis,'BinWidth',50,'Normalization','probability');
xlim([0 600])
set(gca,'FontSize',20,'FontWeight','Bold')
xlabel('Length of principal axis (nm)','FontSize',30)
%ylabel('Frequency','FontSize',30)
set(gcf,'position',[100 100 800 600])

figure;
h2=histogram(short_axis,'BinWidth',50,'Normalization','probability');
xlim([0 600])
set(gca,'FontSize',20,'FontWeight','Bold')
xlabel('Length of auxiliary axis (nm)','FontSize',30)
%ylabel('Frequency','FontSize',30)
set(gcf,'position',[100 100 800 600])

save(strcat(FilePath, '3DNanodomain_', OutputFile, '.mat'));
%% plot individual for testing
% density_thres=0.5;
% index=23;
% figure;
% scatter(Homer_ROI{index}(:,3),Homer_ROI{index}(:,4),'g')
% hold on;
% Output=scatplot(AMPAR_ROI{index}(:,3),AMPAR_ROI{index}(:,4));
% axis equal
% view(-360,-90);
% xlim([center_ROI_new(index,1)-box_size center_ROI_new(index,1)+box_size]);
% ylim([center_ROI_new(index,2)-box_size center_ROI_new(index,2)+box_size]);
% 
% % cluster analysis 
% epsilon=75; % these values are based on Choquet's J neuro paper + empirically obtained. more evidence required.
% minpts=30; % especoa;;u for minpts
% [idx,corepts]=dbscan(AMPAR_ROI{index}(:,3:4),epsilon,minpts);
% data_temp=AMPAR_ROI{index};
% for i=1:length(data_temp);
%     if idx(i)==-1;
%         data_temp(i,:)=zeros(1,12);
%     else
%     end
% end
% 
% data_filtered=data_temp(data_temp(:,9)~=0,:);
% idx_filtered=idx(idx~=-1);
% num_peaks=length(unique(idx_filtered));
% gmmodel=fitgmdist([data_filtered(:,3),data_filtered(:,4)],num_peaks,'CovarianceType','full');
% size_ND=zeros(num_peaks,2);
% for i=1:num_peaks;
%     size_ND(i,1)=2.3*sqrt(gmmodel.Sigma(1,1,i));
%     size_ND(i,2)=2.3*sqrt(gmmodel.Sigma(2,2,i));
% end
% long_axis=max(size_ND')';
% short_axis=min(size_ND')';
% figure;gscatter(data_filtered(:,3),data_filtered(:,4),idx_filtered);axis equal;
% legend('off')
% axis equal
% view(-360,-90);
% xlim([center_ROI_new(index,1)-box_size center_ROI_new(index,1)+box_size]);
% ylim([center_ROI_new(index,2)-box_size center_ROI_new(index,2)+box_size]);
% haxis = gca;
% xlim1 = haxis.XLim;
% ylim1 = haxis.YLim;
% d = (max([xlim1 ylim1])-min([xlim1 ylim]))/1000;
% [X1Grid,X2Grid] = meshgrid(xlim1(1):d:xlim1(2),ylim1(1):d:ylim1(2));
% hold on
% contour(X1Grid,X2Grid,reshape(pdf(gmmodel,[X1Grid(:) X2Grid(:)]),...
%     size(X1Grid,1),size(X1Grid,2)),20);
% 
% % quantification of number
% number_ND=length(unique(idx_filtered));
