
clear;%clear
close all;% Close all figures
% 
path = 'E:\Data\BidentateQD\NewOHsQD\190607_SA-Atto647N_monothiolQD_temperature\3\4\';

if ~exist('fileName1','var')|| isempty(fileName1)
    [userfilein, userdirin]=uigetfile({
        '*.mat','Data file (*.mat)';...
        '*.*','All Files (*.*)'},'Select the trace_R file to process',...
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
        '*.*','All Files (*.*)'},'Select the diffusion_coefficient file to process',...
        path);
    fileName2=fullfile(userdirin,userfilein);
else
    if ~exist(fileName2,'file')
        fprintf('File not found: %s\n',fileName2);
        return;
    end
end
% path = 'E:\experiment\20161202_NMDAR_TRACKING\SA\1nM\FM\1\fPALM\diffusion\';
% 
% load('E:\experiment\20161202_NMDAR_TRACKING\SA\1nM\FM\1\fPALM\diffusion\trace_R.mat')
% load('E:\experiment\20161202_NMDAR_TRACKING\SA\1nM\FM\1\fPALM\diffusion\diffusion_coefficients.mat')
% tr = trace_R;
% s = size(tr);
% s1 = s(1);
% 
% d = Dif_track_W_ID;
% p = size(d);
% p1 = p(1);

% path = 'F:\experiment\NR2A+NR2B_tracking\20180816_DIV15_NR2A_SA_Atto647N\1\test between perisynaptic and synaptic\';
% file1 = 'trace_R.mat';
% file2 = 'diffusion_coefficients_605.mat';
% fname1 = fullfile(path,file1);
% fname2 = fullfile(path,file2);
load(fileName1);
tr = trace_R;
s = size(tr);
s1 = s(1);

load(fileName2);
d = Dif_track_W_ID;
p = size(d);
p1 = p(1);
%%
output= zeros(p1,2);
output2= zeros(p1,4);

for k=1:p1;
   
    output(k,:) = [tr(d(k,2),6) d(k,1)];

    
end

output2(:,1)=output(:,1);
output2(:,2)=output(:,2);

output2(:,3) = log10(output(:,1)/1000);
output2(:,4) = log10(output(:,2));
userfilein_split  = strsplit(userfilein, '.');
matched = strcat(userdirin, userfilein_split(1), '-diffusion_trace.xlsx');
xlswrite(matched{1}, output2);