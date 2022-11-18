clc;clear;
FilePath='E:\Data\DNA-PAINT(Live-cell)\SPT_tracking\cLTP_tracking\210311_DNA-PAINT_cLTP_atRT\Cell2_cLTP_RT_1x\Results';

if ~exist('FileName1','var')|| isempty(FileName1)
    [userfilein, userdirin]=uigetfile({
         '*.txt','Data file (*.txt)';...
        '*.*','All Files (*.*)'},'Select result tracking files to process',...
        FilePath, 'MultiSelect','on');
    FileName1=fullfile(userdirin,userfilein);
else
    if ~exist(FileName1,'file')
        fprintf('File not found: %s\n',FileName1);
        return;
    end
end
Traj=cell(1,length(FileName1));
for i=1:length(FileName1);
    Traj{i}=textread(FileName1{i});
end
Traj_len=cell(1,length(FileName1));
for j=1:length(FileName1);
    [Traj_len{j}(:,1),Traj_len{j}(:,2)]=groupcounts(Traj{j}(:,5));
end

All_length=[];
for k=1:length(FileName1);
    All_length=[All_length; Traj_len{k}(:,1)];
end