% converting ThunderSTORM file to quickPALM format.
clear; close all; clc;
pixel_size=160; 
path = 'E:\Data\BidentateQD\NewOHsQD\190607_SA-Atto647N_monothiolQD_temperature\2\1\TestingThunderSTORM\';
if ~exist('fileName1','var')|| isempty(fileName1)
    [userfilein, userdirin]=uigetfile({
        '*.csv','Data file (*.csv)';...
        '*.*','All Files (*.*)'},'Select the PALM file to process',...
        path);
    fileName1=fullfile(userdirin,userfilein);
else
    if ~exist(fileName1,'file')
        fprintf('File not found: %s\n',fileName1);
        return;
    end
end

tSTORM=csvread([userdirin userfilein],1,0);
qPALM = [tSTORM(:,1) tSTORM(:,8) tSTORM(:,3) tSTORM(:,4) tSTORM(:,3).*pixel_size tSTORM(:,4).*pixel_size  tSTORM(:,5) tSTORM(:,6)./2 tSTORM(:,6)./2 tSTORM(:,7)./2 tSTORM(:,7)./2 tSTORM(:,6) tSTORM(:,7) tSTORM(:,9) tSTORM(:,2)];
xlswrite([path userfilein(1:end-4) '_qPALM.xlsx'], qPALM);