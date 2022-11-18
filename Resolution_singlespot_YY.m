% calculate resolution of PSF.
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
x=cell(1,length(data1));
y=cell(1,length(data1));
del_x=zeros(1,length(data1));
del_y=zeros(1,length(data1));

for i=1:length(data1);
    x{i}=data1{i}(:,3)-mean(data1{i}(:,3));
    y{i}=data1{i}(:,4)-mean(data1{i}(:,4));
    GMModel=fitgmdist([x{i},y{i}],1);
    del_x(i)=sqrt(GMModel.Sigma(1,1));
    del_y(i)=sqrt(GMModel.Sigma(2,2));
end
del_x=del_x';
del_y=del_y';

