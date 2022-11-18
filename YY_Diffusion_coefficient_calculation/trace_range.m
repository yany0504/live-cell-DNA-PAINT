

% This script calculates the range of traces in x, y and z direction
clear;%clear
close all;% Close all figures

%% Read file
if ~exist('fileName1','var')|| isempty(fileName1)
    [userfilein, userdirin]=uigetfile({
        '*.txt','Data file (*.txt)';...
        '*.*','All Files (*.*)'},'Select the result_tracking file to process',...
        'E:\Data\BidentateQD\NewOHsQD\190607_SA-Atto647N_monothiolQD_temperature\3\4\');
    fileName1=fullfile(userdirin,userfilein);
else
    if ~exist(fileName1,'file')
        fprintf('File not found: %s\n',fileName1);
        return;
    end
end
result=textread(fileName1); % Tracking result


%% Group trace according to ID and arrange in cells (unit nm)
x=result(:,1);
y=result(:,2);
z=result(:,3);
frame=result(:,4);
particle=result(:,5);
d=length(particle); %total row number
ptotal=particle(d);%total particle number
for k=1:ptotal
    clear x1 y1 z1 frame1 particle1 dp 
        dp=find(particle==k); %dimension of matrix of particle# K 
        dplength=length(dp);
        dpmax=max(dp);
        dpmin=min(dp);
        m(k)=k; %test k value
        x1=x(dpmin:dpmax);
        y1=y(dpmin:dpmax);
        z1=z(dpmin:dpmax);
        frame1=frame(dpmin:dpmax);
        particle1=particle(dpmin:dpmax);
        trace{k}=[x1,y1,z1,frame1,particle1]; 
        trace_R(k,1)=max(x1)-min(x1); % trace range in X
        trace_R(k,2)=max(y1)-min(y1); % trace range in Y
        trace_R(k,3)=max(z1)-min(z1); % trace range in z
        trace_R(k,4)=max(frame1)-min(frame1);
        trace_R(k,5)=k;    % trace ID       
end
figure
trace_R(:,6)=max(trace_R(:,1:3)');
hist(log10(trace_R(:,6)/1000),48) % trace_r in um
save(strcat(userdirin,'trace_R.mat'),'trace_R'); % trace_R: [X Y Z frame ID max(X Y Z)]
saveas(gcf,strcat(userdirin,'hist_traceRange'),'fig');
clear x y z frame particle k d m x1 y1 z1 frame1 dp ptotal dplength dpmax dpmin
