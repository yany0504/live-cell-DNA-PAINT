%%
clc;clear;close all;
FilePath='E:\Data\DNA-PAINT(Live-cell)\kinetics_calibration\210426_DNA-PAINT_calib_ACSF_1x5x_16xR1_R1_8nt\PSF';

if ~exist('FileName1','var')|| isempty(FileName1)
    [userfilein, userdirin]=uigetfile({
         '*.csv','Data file (*.csv)';...
        '*.*','All Files (*.*)'},'Select single PSF file to process',...
        FilePath, 'MultiSelect','on');
    FileName1=fullfile(userdirin,userfilein);
else
    if ~exist(FileName1,'file')
        fprintf('File not found: %s\n',FileName1);
        return;
    end
end

if iscell(FileName1);
    
    for i=1:length(FileName1);
        data1{i}=xlsread(FileName1{i});
    end
else
    data1=xlsread(FileName1);
end

photons=[];
accuracy=[];
for i=1:length(data1);
    photons=[photons;data1{i}(:,6)];
    accuracy=[accuracy;data1{i}(:,9)];
end
