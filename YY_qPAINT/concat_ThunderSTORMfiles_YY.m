%% concatenating multiple measurement with time step
clear;close all;clc;
path_output = 'H:\Data_from07222021\DNA-PAINT(Live-cell)\SPT_tracking\cLTP_tracking\211009_liveDNA-PAINT_cLTP_noTRF\Cell3_cLTP\';
outputName='Cell3_cLTP';
picasso=0;
exposure=0.05; % second
timestep = 600; % second
dimension=3;

interval=timestep/exposure;

%%
[userfilein, userdirin]=uigetfile({
    '*.csv*','Data file (*.csv)';...
    '*.*','All Files (*.*)'},'Select the ThunderSTORM files to process',...
    path_output,'MultiSelect','on');
fileName_TS=fullfile(userdirin,userfilein);
if iscell(fileName_TS);
    for i=1:length(fileName_TS);
        try
            ThunderSTORM{i,1}=dlmread(fileName_TS{i},'\t',1,0);
        catch exception
            ThunderSTORM{i,1}=xlsread(fileName_TS{i});
        end
        if picasso ==1
            ThunderSTORM{i,1}(:,2)=ThunderSTORM(:,2)+1;
        else
        end
    end
    numberImages=length(fileName_TS);
    prompt_frame=cell(2*numberImages,1);
else
    try
        ThunderSTORM1 = dlmread(fileName_TS,'\t',1,0);
    catch exception
        ThunderSTORM1 = xlsread(fileName_TS);
    end
    
    if picasso ==1
        ThunderSTORM1(:,2)=ThunderSTORM1(:,2)+1;
    else
    end
    numberImages=1;
    prompt_frame=cell(2*numberImages,1);
    ThunderSTORM=cell(1,1);
    ThunderSTORM{1,1}=ThunderSTORM1;
    fileName_TS1=fileName_TS;
    fileName_TS=cell(1,1);
    fileName_TS{1,1}=fileName_TS1;
end
%% adding time steps and concat
for i=1:length(ThunderSTORM);
    ThunderSTORM{i}(:,2)=ThunderSTORM{i}(:,2)+(i-1)*interval;
    
end
ThunderSTORM_concat=ThunderSTORM{1};
for j=2:length(ThunderSTORM);
    ThunderSTORM_concat=[ThunderSTORM_concat;ThunderSTORM{j}];
end
%% save concatenated file.

headers_3D={'id','frame','x [nm]','y [nm]','z [nm]','sigma1 [nm]','sigma2 [nm]','intensity [photon]','offset [photon]','bkgstd [photon]','uncertainty_xy [nm]','uncertainty_z [nm]'};
headers_2D={'id','frame','x [nm]','y [nm]','sigma [nm]','intensity [photon]','offset [photon]','bkgstd [photon]','uncertainty_xy [nm]'};

if dimension == 3;
    T=array2table(ThunderSTORM_concat,'VariableNames',headers_3D);
    writetable(T,strcat([path_output outputName '_concat.csv']));
elseif dimension == 2;
    
    T=array2table(ThunderSTORM_concat,'VariableNames',headers_2D);
    writetable(T,strcat([path_output outputName '_concat.csv']));
else
    fprintf('Dimension is wrong');
end


