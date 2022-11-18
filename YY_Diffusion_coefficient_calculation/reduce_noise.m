function reduce_noise(fileName1,fileName2, fileName3)
clc;%clear
close all;

path='F:\DATA_3\cell imaging\NMDA and AMPA receptors\AMPAR\20151027_AMPAR_bigqdot(625nm)\8'
if ~exist('fileName1','var')|| isempty(fileName1)
    [userfilein, userdirin]=uigetfile({
        '*.xlsx','Data file (*.xlsx)';...
        '*.*','All Files (*.*)'},'Select the synaptotagmin file to process',...
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
        '*.xlsx','Data file (*.xlsx)';...
        '*.*','All Files (*.*)'},'Select the receptor file to process',...
        path);
    fileName2=fullfile(userdirin,userfilein);
else
    if ~exist(fileName2,'file')
        fprintf('File not found: %s\n',fileName2);
        return;
    end
end
if ~exist('fileName3','var')|| isempty(fileName3)
    [userfilein, userdirin]=uigetfile({
        '*.xlsx','Data file (*.xlsx)';...
        '*.*','All Files (*.*)'},'Select the homer1 file to process',...
        path);
    fileName3=fullfile(userdirin,userfilein);
else
    if ~exist(fileName3,'file')
        fprintf('File not found: %s\n',fileName3);
        return;
    end
end
%% synapse_2.m for Green
data=xlsread(fileName1);
palmX=data(:,5);
palmY=data(:,6);
palmZ=data(:,7);
data_Syn=[palmX, palmY, palmZ];

% Finding clusters
[Class,type]=dbscan_conservative(data_Syn,50,500); 
% Make new matrix  
Syn=[palmX,palmY,palmZ,Class',type'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seperate and plot clusters
x=Syn(:,1);
y=Syn(:,2);
z=0.79*Syn(:,3); 
particle=Syn(:,4);
PtType=Syn(:,5);
d=length(particle); %total row number
ptotal=max(particle);%total cluster number
numArrays=ptotal;
Synapses1=cell(numArrays,1);
figure
for k=1:ptotal
        dp=find(particle==k); % matrix indcis of points in cluster #k 
        dplength=length(dp);
        %dpmax=max(dp);
        %dpmin=min(dp);
        m(k)=k; %test k value
        x1=x(dp);
        y1=y(dp);
        z1=z(dp);
        PtType1=PtType(dp);
        Syn1=[x1,y1,z1,PtType1];
        Synapses1{k}=Syn1;
        %eval(['Synapse_' num2str(k) '=[x1,y1,z1,PtType1]']);% creating submatrix for particle# k sub_k 
        scatter3(x1,y1,z1,10,'.');
        hold all
        clear dp dplength dpmax dpmin x1 y1 z1 PtType1
end

hold off
clear x y z PtType particle k d m x1 y1 z1 dp Syn 
axis equal
% load('D:\mydocument\MATLAB\for codes\Synapses2.mat')
% load('D:\mydocument\MATLAB\for codes\Synapses1.mat')
% userdirin='D:\mydocument\MATLAB\for codes\';
% ptotal=size(Synapses1,1);
fid=fopen([userdirin 'Synap_location.txt'],'w');
for k=1:ptotal
    for m=1:size(Synapses1{k},1)
    fprintf(fid,'%f %f %f %d\n',Synapses1{k}(m,:));
    end
end
Synap = Synapses1;
fclose(fid);
save(strcat(userdirin,'synap.mat'),'Synap');
%% synapse_2.m for red
data=xlsread(fileName2);

palmX=data(:,5);
palmY=data(:,6);
palmZ=data(:,7);
data_Syn=[palmX, palmY, palmZ];

% Finding clusters
[Class,type]=dbscan_conservative(data_Syn,50,500); 
% Make new matrix  
Syn=[palmX,palmY,palmZ,Class',type'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seperate and plot clusters
x=Syn(:,1);
y=Syn(:,2);
z=0.79*Syn(:,3); 
particle=Syn(:,4);
PtType=Syn(:,5);
d=length(particle); %total row number
ptotal=max(particle);%total cluster number
numArrays=ptotal;
Synapses2=cell(numArrays,1);
figure
for k=1:ptotal
        dp=find(particle==k); % matrix indcis of points in cluster #k 
        dplength=length(dp);
        %dpmax=max(dp);
        %dpmin=min(dp);
        m(k)=k; %test k value
        x1=x(dp);
        y1=y(dp);
        z1=z(dp);
        PtType1=PtType(dp);
        Syn1=[x1,y1,z1,PtType1];
        Synapses2{k}=Syn1;
        %eval(['Synapse_' num2str(k) '=[x1,y1,z1,PtType1]']);% creating submatrix for particle# k sub_k 
        scatter3(x1,y1,z1,10,'.');
        hold all
        clear dp dplength dpmax dpmin x1 y1 z1 PtType1
end

hold off
clear x y z PtType particle k d m x1 y1 z1 dp Syn 
axis equal
% load('D:\mydocument\MATLAB\for codes\Synapses2.mat')
% load('D:\mydocument\MATLAB\for codes\Synapses1.mat')
% userdirin='D:\mydocument\MATLAB\for codes\';
% ptotal=size(Synapses2,1);
fid=fopen([userdirin 'NMDAR_locations.txt'],'w');
for k=1:ptotal
    for m=1:size(Synapses2{k},1)
    fprintf(fid,'%f %f %f %d\n',Synapses2{k}(m,:));
    end
end

NMDAR = Synapses2;
fclose(fid);
save(strcat(userdirin,'NMDAR.mat'),'NMDAR');

%% synapse_3.m for red
data=xlsread(fileName3);

palmX=data(:,5);
palmY=data(:,6);
palmZ=data(:,7);
data_Syn=[palmX, palmY, palmZ];

% Finding clusters
[Class,type]=dbscan_conservative(data_Syn,50,500); 
% Make new matrix  
Syn=[palmX,palmY,palmZ,Class',type'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seperate and plot clusters
x=Syn(:,1);
y=Syn(:,2);
z=0.79*Syn(:,3); 
particle=Syn(:,4);
PtType=Syn(:,5);
d=length(particle); %total row number
ptotal=max(particle);%total cluster number
numArrays=ptotal;
Synapses3=cell(numArrays,1);
figure
for k=1:ptotal
        dp=find(particle==k); % matrix indcis of points in cluster #k 
        dplength=length(dp);
        %dpmax=max(dp);
        %dpmin=min(dp);
        m(k)=k; %test k value
        x1=x(dp);
        y1=y(dp);
        z1=z(dp);
        PtType1=PtType(dp);
        Syn1=[x1,y1,z1,PtType1];
        Synapses3{k}=Syn1;
        %eval(['Synapse_' num2str(k) '=[x1,y1,z1,PtType1]']);% creating submatrix for particle# k sub_k 
        scatter3(x1,y1,z1,10,'.');
        hold all
        clear dp dplength dpmax dpmin x1 y1 z1 PtType1
end

hold off
clear x y z PtType particle k d m x1 y1 z1 dp Syn 
axis equal
% load('D:\mydocument\MATLAB\for codes\Synapses2.mat')
% load('D:\mydocument\MATLAB\for codes\Synapses1.mat')
% userdirin='D:\mydocument\MATLAB\for codes\';
% ptotal=size(Synapses2,1);
fid=fopen([userdirin 'homer1_locations.txt'],'w');
for k=1:ptotal
    for m=1:size(Synapses3{k},1)
    fprintf(fid,'%f %f %f %d\n',Synapses3{k}(m,:));
    end
end

Homer = Synapses3;
fclose(fid);
save(strcat(userdirin,'Homer.mat'),'Homer');

