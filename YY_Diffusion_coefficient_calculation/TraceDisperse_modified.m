clear all;%clear
close all;% Close all figures

path = 'E:\Data\BidentateQD\NewOHsQD\190607_SA-Atto647N_monothiolQD_temperature\3\4\';
if ~exist('fileName1','var')|| isempty(fileName1)
    [userfilein, userdirin]=uigetfile({
        '*.txt','Data file (*.txt)';...
        '*.*','All Files (*.*)'},'Select the result_tracking to process',...
        path);
    fileName1=fullfile(userdirin,userfilein);
else
    if ~exist(fileName1,'file')
        fprintf('File not found: %s\n',fileName1);
        return;
    end
end
if ~exist('fileName2','var')|| isempty(fileName2)
    [userfilein, userdirin]=uigetfile({
        '*.mat','Data file (*.mat)';...
        '*.*','All Files (*.*)'},'Select the Syanpse1 file to process',...
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
result_tracking=textread(fileName1);
S_Synapses1=load(fileName2);
Synapses1=S_Synapses1.Synapses1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step2
%%%%%%%%%%%%%
%%%%%%trace_center
% This script finds the center of traces abtained from tracking.m
%===============================
% load data from tracking.m
x=result_tracking(:,1);
y=result_tracking(:,2);
z=result_tracking(:,3);
frame=result_tracking(:,4);
particle=result_tracking(:,5);
d=length(particle); %total row number
ptotal=particle(d); %total particle number
for n=1:ptotal
    dp=find(particle==n); %dimension of matrix of particle# K
%     dplength=length(dp);
    dpmax=max(dp);
    dpmin=min(dp);
    m(n)=n; %test k value
    x1=x(dpmin:dpmax);
    y1=y(dpmin:dpmax);
    z1=z(dpmin:dpmax);
    pt1=particle(dpmin:dpmax);
    %trace=[x1,y1,z1];
    frame1=frame(dpmin:dpmax);
    traces{n}=[x1,y1,z1,frame1,pt1];% creating a cell array for particle# k traces{k}
    traces_center(n,1)=mean(x1);
    traces_center(n,2)=mean(y1);
    traces_center(n,3)=mean(z1);
    

% save(strcat(userdirin,'traces.mat'),'traces');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  nSyn=length(Synapses1);
  
  % syn_dis=2000;  % unit in nm
  
 
    
      for rg=1:nSyn
        syn_rg=Synapses1{rg}; % load data
        syn_length=length(syn_rg);
        % Find center of the synapse
        clear trace_nearby syn_traces
        %radii=1; % Edit by Anne 11/26/2013 
        Ctr0=mean(syn_rg);%ORIGINAL         [Ctr00,Sig00]=subclust(syn_rg,radii); % Find center of a subcluster and the range of influence of the center
        %%make sure Ctr0 has only 1 row
        dist=zeros(length(traces{1}));
        % Find traces close to the synapse
        syn_traces(rg,1)=(Ctr0(1)-traces_center(n,1)).*(Ctr0(1)-traces_center(n,1));
        syn_traces(rg,2)=(Ctr0(2)-traces_center(n,2)).*(Ctr0(2)-traces_center(n,2));
        syn_traces(rg,3)=(Ctr0(3)-traces_center(n,3)).*(Ctr0(3)-traces_center(n,3));
        distance(rg)=sqrt(syn_traces(rg,1)+syn_traces(rg,2)+syn_traces(rg,3))/1000;
        dist = [distance];
dist2 = dist';
mini=min(dist2);
        t= find(dist2 == mini);
        a=Synapses1{t};
      syn_cent0= mean(a);
     Syn_ctrX = syn_cent0(1);
     Syn_ctrY = syn_cent0(2);
     Syn_ctrZ = syn_cent0(3);
       
       
        TS_X=(traces{n}(:,1)-Syn_ctrX)/1000; % Trace to synapse distance delta X
        TS_Y=(traces{n}(:,2)-Syn_ctrY)/1000; 
        TS_Z=(traces{n}(:,3)-Syn_ctrZ)/1000;
      
      
      TS_dist{n}=sqrt(TS_X.*TS_X+TS_Y.*TS_Y+TS_Z.*TS_Z);
        %syn_nearby{rg}=find(syn_traces(:,4)<syn_dis);
        %syn_nearby_t{rg}=find(syn_traces(:,4)<syn_dis)'; % syn_nearby_t is for coverting into traceID array
        %Syn_ctr(rg,:)=Ctr0(1:3);% Add by Anne 11/26/2013
      end
    
   % again synapess center

      
end

 save(strcat(userdirin,'TS_dist.mat'),'TS_dist');
figure; hist(cell2mat(TS_dist'),48) 
saveas(gcf,strcat(userdirin,'TS_dist'),'fig');

distance=TS_dist';

D=cell2mat(distance);

% L=size(D_table);



% fclose(fid);
xlswrite(strcat(path,'TS_dist_cal.xlsx'), D);

%     save(strcat(userdirin,'Syn_ctr.mat'),'Syn_ctr');
%     syn_nearby0=find(~cellfun(@isempty,syn_nearby)); % Synapses ID with traces around
  
  

