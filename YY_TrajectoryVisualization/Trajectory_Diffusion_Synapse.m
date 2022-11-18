% Plot trajectories (function of diffusion coefficient) and synapses.
%% Step 1. loading files
clc;clear;close all;
FilePath = 'H:\Data_from07222021\DNA-PAINT(Live-cell)\SPT_tracking\cLTP_tracking\211007_liveDNA-PAINT_cLTP\Cell2_cLTP\';
OutputFile='test';
if ~exist('FileName_resulttracking','var')|| isempty(FileName_resulttracking)
    [userfilein, userdirin]=uigetfile({
         '*.txt','Data file (*.txt)';...
        '*.*','All Files (*.*)'},'Select result_tracking file to process',...
        FilePath, 'MultiSelect','on');
    FileName_resulttracking=fullfile(userdirin,userfilein);
else
    if ~exist(FileName_resulttracking,'file')
        fprintf('File not found: %s\n',FileName_resulttracking);
        return;
    end
end
ResultTracking=dlmread(FileName_resulttracking);
if ~exist('FileName_Synapse','var')|| isempty(FileName_Synapse)
    [userfilein, userdirin]=uigetfile({
         '*.mat','Data file (*.mat)';...
        '*.*','All Files (*.*)'},'Select Synapses1 file to process',...
        FilePath, 'MultiSelect','on');
    FileName_Synapse=fullfile(userdirin,userfilein);
else
    if ~exist(FileName_Synapse,'file')
        fprintf('File not found: %s\n',FileName_Synapse);
        return;
    end
end
load(FileName_Synapse);
if ~exist('FileName_Diff_W_ID','var')|| isempty(FileName_Diff_W_ID)
    [userfilein, userdirin]=uigetfile({
         '*.mat','Data file (*.mat)';...
        '*.*','All Files (*.*)'},'Select diffusion_coefficients_xxx.mat file to process',...
        FilePath, 'MultiSelect','on');
    FileName_Diff_W_ID=fullfile(userdirin,userfilein);
else
    if ~exist(FileName_Diff_W_ID,'file')
        fprintf('File not found: %s\n',FileName_Diff_W_ID);
        return;
    end
end

load(FileName_Diff_W_ID);
Dif_track_W_ID(:,1)=log10(Dif_track_W_ID(:,1));
ID=Dif_track_W_ID(:,2);

%% Step 2. generate colormap for diffusion
Diffmin=-5;
Diffmax=-0.5;
Diffstep=11;
figure;
ax = axes;
cmap=colormap(jet(Diffstep));
c=colorbar('Ticks',[0,1],'TickLabels',{['10^' '{' num2str(Diffmin) '}'], ['10^' '{' num2str(Diffmax) '}']},'FontSize',20,'Location','west');
c.Label.String='Diffusion coefficient (\mum^2/s)';
ax.Visible = 'off';
Diff_range=linspace(Diffmin+abs(Diffmin+Diffmax)/(2*Diffstep),Diffmax-abs(Diffmin+Diffmax)/(2*Diffstep),Diffstep);

%% Step 3. Determine color code. find closest color code for diffusion coefficients
ColorCode=zeros(length(Dif_track_W_ID),1);
for i=1:length(Dif_track_W_ID)
    [val, ColorCode(i)]=min(abs(Diff_range-Dif_track_W_ID(i,1)));
end
%% Step 4. Find the ID of trace_R in Dif_track_W_ID(:,2) and assign corresponding colorcode
for i=1:length(ResultTracking)
    if ismember(ResultTracking(i,5),ID)==0
        ResultTracking(i,:)=zeros(1,5);
    else
    end  
end
ResultTrackingFiltered=ResultTracking(any(ResultTracking,2),:);
ColorID=zeros(length(ResultTrackingFiltered),1);
for i=1:length(ResultTrackingFiltered)
    ColorID(i)=ColorCode(find(ID==ResultTrackingFiltered(i,5)));
end
x=ResultTrackingFiltered(:,1);
y=ResultTrackingFiltered(:,2);
z=ResultTrackingFiltered(:,3);
t=ResultTrackingFiltered(:,4);
particleID=ResultTrackingFiltered(:,5);


%% Step 5.1. Plot trajectory with synapse ID-labeled side-by-side

traj_diff=[particleID t x y z ColorID];
%traj_diff=[x y z particleID ColorID];
[Uid,~,ix]=unique(traj_diff(:,1));
tally = accumarray(ix, ones(size(ix)));
trajectories = mat2cell(traj_diff, tally, [6]);
figure;
ax1=subplot(1,2,1);
hold on;
% for 2D
for k = 1:size(trajectories,1)
    plot(trajectories{k}(:,3), trajectories{k}(:,4),'Color',cmap(trajectories{k}(1,6),:),'LineWidth',2);
%     h=text(trajectories{k}(1,3),trajectories{k}(1,4),trajectories{k}(1,5),num2str(trajectories{k}(1,1)),'Color','white','FontSize',15);
%     set(h, 'Clipping','on');
end
% for 3D
% for k = 1:size(trajectories,1)
%     plot3(trajectories{k}(:,3), trajectories{k}(:,4),trajectories{k}(:,5),'Color',cmap(trajectories{k}(1,6),:),'LineWidth',2);
% %     h=text(trajectories{k}(1,3),trajectories{k}(1,4),trajectories{k}(1,5),num2str(trajectories{k}(1,1)),'Color','white','FontSize',15);
% %     set(h, 'Clipping','on');
% end

data=Synapses1{1};
for i=2:length(Synapses1)
    data=cat(1,data,Synapses1{i});
end
%densityScatter([data(:,1),data(:,2),data(:,3)]','colMap',bone(10),'filterFactor',0.9,'markerSize',50);
scatter(data(:,1),data(:,2),'w','Marker','o');
ax1 = gca;
set(ax1,'Color','k');
axis equal

view(-360,-90);
grid on
% zlim([-600 600])
% ax.GridColor = [0,0,0]; %black grid
ax1.GridColor = [1, 1, 1]; % white grid
% ax1.Position = [0.1 0.1 0.8 0.8];
ax1.YAxis.FontSize=30;
ax1.XAxis.FontSize=30;
ax1.ZAxis.FontSize=30;
set(gca,'FontWeight','bold');
xlabel('x (nm)','FontSize',40,'fontweight','bold');
ylabel('y (nm)','FontSize',40,'fontweight','bold');
zlabel('z (nm)','FontSize',40,'fontweight','bold');

ax2=subplot(1,2,2);
hold on;

% for k = 1:size(trajectories,1)
%     plot3(trajectories{k}(:,3), trajectories{k}(:,4),trajectories{k}(:,5),'Color',cmap(trajectories{k}(1,6),:),'LineWidth',2);
%     h=text(trajectories{k}(1,3),trajectories{k}(1,4),trajectories{k}(1,5),num2str(trajectories{k}(1,1)),'Color','white','FontSize',15);
%     set(h, 'Clipping','on');
% end
% for 2D
for k = 1:size(trajectories,1)
    plot(trajectories{k}(:,3), trajectories{k}(:,4),'Color',cmap(trajectories{k}(1,6),:),'LineWidth',2);
    h=text(trajectories{k}(1,3),trajectories{k}(1,4),num2str(trajectories{k}(1,1)),'Color','white','FontSize',15);
    set(h, 'Clipping','on');
end
hold off;
linkaxes([ax1 ax2], 'xy');
ax2 = gca;
set(ax2,'Color','k');
axis equal
view(-360,-90);
grid on
% ax.GridColor = [0,0,0]; %black grid
ax2.GridColor = [1, 1, 1]; % white grid
% ax2.Position = [0.1 0.1 0.8 0.8];
ax2.YAxis.FontSize=30;
ax2.XAxis.FontSize=30;
ax2.ZAxis.FontSize=30;
set(gca,'FontWeight','bold');
xlabel('x (nm)','FontSize',40,'fontweight','bold');
ylabel('y (nm)','FontSize',40,'fontweight','bold');
zlabel('z (nm)','FontSize',40,'fontweight','bold');
xlim([32000 33500]);
ylim([26000 27500]);
%% Step 5.2. Plot trajectory

traj_diff=[particleID t x y z ColorID];
%traj_diff=[x y z particleID ColorID];
[Uid,~,ix]=unique(traj_diff(:,1));
tally = accumarray(ix, ones(size(ix)));
trajectories = mat2cell(traj_diff, tally, [6]);
figure;
hold on;
% % for 2D
% for k = 1:size(trajectories,1)
%     plot(trajectories{k}(:,3), trajectories{k}(:,4),'Color',cmap(trajectories{k}(1,6),:),'LineWidth',2);
% %     h=text(trajectories{k}(1,3),trajectories{k}(1,4),trajectories{k}(1,5),num2str(trajectories{k}(1,1)),'Color','white','FontSize',15);
% %     set(h, 'Clipping','on');
% end
% for 3D
for k = 1:size(trajectories,1)
    plot3(trajectories{k}(:,3), trajectories{k}(:,4),trajectories{k}(:,5),'Color',cmap(trajectories{k}(1,6),:),'LineWidth',2);
%     h=text(trajectories{k}(1,3),trajectories{k}(1,4),trajectories{k}(1,5),num2str(trajectories{k}(1,1)),'Color','white','FontSize',15);
%     set(h, 'Clipping','on');
end

data=Synapses1{1};
for i=2:length(Synapses1)
    data=cat(1,data,Synapses1{i});
end
%densityScatter([data(:,1),data(:,2),data(:,3)]','colMap',bone(10),'filterFactor',0.9,'markerSize',50);
scatter3(data(:,1),data(:,2),data(:,3),'w','Marker','.');
ax1 = gca;
set(ax1,'Color','k');
axis equal

view(-360,-90);
grid on
% zlim([-600 600])
% ax.GridColor = [0,0,0]; %black grid
ax1.GridColor = [1, 1, 1]; % white grid
% ax1.Position = [0.1 0.1 0.8 0.8];
ax1.YAxis.FontSize=30;
ax1.XAxis.FontSize=30;
ax1.ZAxis.FontSize=30;
set(gca,'FontWeight','bold');
xlabel('x (nm)','FontSize',40,'fontweight','bold');
ylabel('y (nm)','FontSize',40,'fontweight','bold');
zlabel('z (nm)','FontSize',40,'fontweight','bold');
zlim([-600 600])

%% filter out trajectories diffcoeff<10^-2
diff_thres=-1.95;
trajectories_slow=cell(length(trajectories),1);
for i=1:length(trajectories);
    if Dif_track_W_ID(i,1)<=diff_thres;
        trajectories_slow{i,1}=trajectories{i,1};
    else
    end
end
Dif_slow_W_ID=Dif_track_W_ID(Dif_track_W_ID(:,1)<=diff_thres,:);
trajectories_slow=trajectories_slow(~cellfun('isempty',trajectories_slow));

figure;
hold on;
% for 2D
for k = 1:size(trajectories_slow,1)
    plot(trajectories_slow{k}(:,3), trajectories_slow{k}(:,4),'Color',cmap(trajectories_slow{k}(1,6),:),'LineWidth',2);
end
data=Synapses1{1};
for i=2:length(Synapses1)
    data=cat(1,data,Synapses1{i});
end
%densityScatter([data(:,1),data(:,2),data(:,3)]','colMap',bone(10),'filterFactor',0.9,'markerSize',50);
scatter(data(:,1),data(:,2),'w','Marker','o');
ax1 = gca;
set(ax1,'Color','k');
axis equal

view(-360,-90);
grid on
% zlim([-600 600])
% ax.GridColor = [0,0,0]; %black grid
ax1.GridColor = [1, 1, 1]; % white grid
% ax1.Position = [0.1 0.1 0.8 0.8];
ax1.YAxis.FontSize=30;
ax1.XAxis.FontSize=30;
ax1.ZAxis.FontSize=30;
set(gca,'FontWeight','bold');
xlabel('x (nm)','FontSize',40,'fontweight','bold');
ylabel('y (nm)','FontSize',40,'fontweight','bold');
zlabel('z (nm)','FontSize',40,'fontweight','bold');

%% setting visualization
xticks('auto');
yticks('auto');


% x0=100;
% y0=100;
% width=1000;
% height=600;
% set(gcf,'position',[x0,y0,width,height])
% view(-57.6283, 11.7389);
% ylim([17400 18000]);
% xlim([20200 20700]); 
% zlim([-350 -50]);
% ax.YAxis.Exponent = 0;
% set(gca,'Xtick',0:200:80000);
% set(gca,'Ytick',0:200:80000);
% xtickangle(-17)
% ax.XAxis.Exponent = 0;
% 
% ytickangle(17)
% ax.XLabel.Position = [1.99e+04 1.7e+04 -300];
% ax.XLabel.Rotation =22;
% ax.YLabel.Position = [1.96e+04 1.74e+04 -300];
% ax.YLabel.Rotation =-16;

%% preparation for animated figure
traj_num=2;
index=find(Uid==traj_num);
range=5000; %nm

data1=data(data(:,1)>=(mean(trajectories{index}(:,3))-range)&data(:,1)<(mean(trajectories{index}(:,3))+range)&...
    data(:,2)>=(mean(trajectories{index}(:,4))-range)&data(:,2)<(mean(trajectories{index}(:,4))+range),:);
h = figure;
densityScatter([data1(:,1),data1(:,2),data1(:,3)]','colMap',bone(10),'filterFactor',0.9,'markerSize',100);
hold on


ax = gca;
set(ax,'Color','k');
axis equal
grid on
% ax.GridColor = [0,0,0]; %black grid
ax.GridColor = [1, 1, 1]; % white grid
ax.Position = [0.1 0.1 0.8 0.8];
ax.YAxis.FontSize=30;
ax.XAxis.FontSize=30;
ax.ZAxis.FontSize=30;
set(gca,'FontWeight','bold');
xlabel('x (nm)','FontSize',40,'fontweight','bold');
ylabel('y (nm)','FontSize',40,'fontweight','bold');
zlabel('z (nm)','FontSize',40,'fontweight','bold');

xlim([mean(trajectories{index}(:,3))-range mean(trajectories{index}(:,3))+range]);
ylim([mean(trajectories{index}(:,4))-range mean(trajectories{index}(:,4))+range]);
zlim([-300 300]);
view(3);
%% animated figure gogogo!
% if you need...
xlim([33500 35500]);
ylim([11000 15000]);
ax.Position = [0.11 0.1 0.8 0.8];
traj_num=1691;
index=find(Uid==traj_num);

for n = 1:size(trajectories{index},1)
    h1=plot3(trajectories{index}(1:n,3), trajectories{index}(1:n,4),trajectories{index}(1:n,5),'Color',cmap(trajectories{index}(1,6),:),'LineWidth',2); 
    h2=scatter3(trajectories{index}(n,3), trajectories{index}(n,4),trajectories{index}(n,5),75,cmap(trajectories{index}(1,6),:),'filled','Linewidth',2);
    drawnow
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if n == 1 
          imwrite(imind,cm,[FilePath OutputFile '_' num2str(traj_num) '_animated.gif'],'gif','DelayTime',0.05, 'Loopcount',inf); 
      else 
          imwrite(imind,cm,[FilePath OutputFile '_' num2str(traj_num) '_animated.gif'],'gif','DelayTime',0.05,'WriteMode','append'); 
      end 
     delete(h1);
     delete(h2);
end

% x0=100;
% y0=100;
% width=1000;
% height=600;
% set(gcf,'position',[x0,y0,width,height])
% view(-57.6283, 11.7389);
% ylim([17400 18000]);
% xlim([20200 20700]); 
% zlim([-350 -50]);
% ax.YAxis.Exponent = 0;
% set(gca,'Xtick',0:200:80000);
% set(gca,'Ytick',0:200:80000);
% xtickangle(-17)
% ax.XAxis.Exponent = 0;
% 
% ytickangle(17)
% ax.XLabel.Position = [1.99e+04 1.7e+04 -300];
% ax.XLabel.Rotation =22;
% ax.YLabel.Position = [1.96e+04 1.74e+04 -300];
% ax.YLabel.Rotation =-16;


%%

xlim([35000 38000]);
ylim([14000 17000]);
zlim([-300 300]);
ax.Position = [0.15 0.1 0.8 0.8];
traj_num=[2287, 2219, 302, 2341, 2338, 1667, 2093];


for i=1:length(traj_num)
    index=find(Uid==traj_num(i));
    n=size(trajectories{index},1);
    h1=plot3(trajectories{index}(1:n,3), trajectories{index}(1:n,4),trajectories{index}(1:n,5),'Color',cmap(trajectories{index}(1,6),:),'LineWidth',2); 
end

saveas(gcf,[FilePath OutputFile '_traj.fig'],'fig');



