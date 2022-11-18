% diffusion coefficient of synaptic (<500 nm) vs. extra
% synaptic (>500 nm)AMPARs
clear; clc; %close all;

%% load distance vs. diffusion file column 1: diffusion coefficient, 2: log(diffusion coefficient) 3: distance 4: ID

FilePath = 'H:\Data_from07222021\DNA-PAINT(Live-cell)\SPT_tracking\Basal_tracking\SPT_results_0701_0708_0715\';
OutputFile='Dist_Diff';
synapse_size=0.5; %micron
if ~exist('FileName1','var')|| isempty(FileName1)
    [userfilein, userdirin]=uigetfile({
         '*.xlsx','Data file (*.xlsx)';...
        '*.*','All Files (*.*)'},'Select the diffusion vs trajectory file to process',...
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

if iscell(data1);
    DiffDist=data1{1};
    for i=2:length(FileName1);
        DiffDist=cat(1,DiffDist,data1{i});
    end
else
    DiffDist=data1;
end
%% Distinguishing synaptic vs. Extrasynaptic
SynapticDiff=[];
ExtrasynapticDiff=[];

for i=1:length(DiffDist);
    if DiffDist(i,3) <= synapse_size;
        SynapticDiff=[SynapticDiff; DiffDist(i,:)];
    elseif DiffDist(i,3) > synapse_size;
        ExtrasynapticDiff=[ExtrasynapticDiff; DiffDist(i,:)];
    end
end

%% parameters
Diffbin=0.1;
Diffmin=-5.0;
Diffmax=0;
FigSize=[100 100 1800 600];
figure;
subplot(1,2,1);
hist_synapticdiff=histogram(SynapticDiff(:,2),'BinWidth',Diffbin,'BinLimits',[Diffmin-0.5*Diffbin,Diffmax+0.5*Diffbin],'Normalization','probability');
xlim([Diffmin-0.5*Diffbin Diffmax+0.5*Diffbin])
title('Synaptic AMPAR','FontSize',25);
ax=gca;
ax.YAxis.FontSize=20;
ax.XAxis.FontSize=20;
xlabel('Diffusion coefficient (\mum^2/s, log)','FontSize',25);
ylabel('Frequency','FontSize',25);


subplot(1,2,2);
hist_extrasynapticdiff=histogram(ExtrasynapticDiff(:,2),'BinWidth',Diffbin,'BinLimits',[Diffmin-0.5*Diffbin,Diffmax+0.5*Diffbin],'Normalization','probability');
xlim([Diffmin-0.5*Diffbin Diffmax+0.5*Diffbin])
title('Extrasynaptic AMPAR','FontSize',25);
ax=gca;
ax.YAxis.FontSize=20;
ax.XAxis.FontSize=20;
xlabel('Diffusion coefficient (\mum^2/s, log)','FontSize',25);
ylabel('Frequency','FontSize',25);

set(gcf, 'Position', FigSize);
saveas(gcf,[FilePath OutputFile 'syn_extrasyn_diffcoeff'],'fig');
