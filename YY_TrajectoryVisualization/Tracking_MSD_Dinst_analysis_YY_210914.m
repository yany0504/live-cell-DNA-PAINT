clear;clc;
tracking=textread('H:\Data_from07222021\DNA-PAINT(Live-cell)\SPT_tracking\Basal_tracking\210708_liveDNA-PAINT_5xR1_16xR1\OnlyForSPT\SPT_results_210708\result_tracking_9_mem3_good10_maxdisp500_k50_eps300.txt');
exposure=0.1;
%%
ind_MSD=cell(max(tracking(:,5)),1);
ind_time=cell(max(tracking(:,5)),1);
Diff=zeros(max(tracking(:,5)),1);
range=zeros(max(tracking(:,5)),1);
instDiff=[];
%%
figure;
for i=1:max(tracking(:,5));
    x=tracking(tracking(:,5)==i,1);
    y=tracking(tracking(:,5)==i,2);
    z=tracking(tracking(:,5)==i,3);
    t=tracking(tracking(:,5)==i,4)*exposure;
    x0=x-x(1);
    y0=y-y(1);
    z0=z-z(1);
    t0=t-t(1);
    subplot(1,3,1);hold on;
    plot3(x0,y0,z0);
    ind_MSD{i}=(x0.^2+y0.^2+z0.^2);
    ind_time{i}=t0;
%     subplot(1,3,2);hold on;
%     plot(t0,ind_MSD{i});
    MSD_sm=smooth(t0,ind_MSD{i},0.2);
    subplot(1,3,2);hold on;
    plot(t0,MSD_sm);
    deriv_msd_s=diff(MSD_sm)/exposure;
    D_inst=abs(deriv_msd_s/(6*10^6));
    subplot(1,3,3);hold on;
    plot(t0(1:end-1)+0.5*exposure,D_inst);
    Diff(i)=mean(D_inst);
    instDiff=[instDiff;D_inst];
    range(i)=sqrt(max(MSD_sm));
end
Diff_log=log10(Diff);
Traj_log=log10(range/1000);
figure;histogram(Diff_log);
figure;histogram(log10(instDiff));
figure;histogram(log10(range/1000));
figure;scatter(log10(Diff),log10(range/1000));
%% average MSD
time_ave=(0:exposure:50)';
msd_tot=zeros(length(time_ave),length(ind_MSD));
for j=1:length(ind_MSD);
    
    
end
%% testing

levelStep=1;
[GMModel_tot,GMModel_diff_tot,GMModel_traj_tot]=Make2Dcontour_MLE(Diff_log,Traj_log,'levelstep',levelStep,'title','Total','PlotStyle','scatter','LineColor','none','colorcode','parula','text','on','histograms','off','Numpopulation',2);
% i=1334
% x=tracking(tracking(:,5)==i,1);
% y=tracking(tracking(:,5)==i,2);
% z=tracking(tracking(:,5)==i,3);
% t=tracking(tracking(:,5)==i,4)*exposure;
% x0=x-x(1);
% y0=y-y(1);
% z0=z-z(1);
% t0=t-t(1);
% subplot(1,3,1);hold on;
% plot3(x0,y0,z0);
% ind_MSD{i}=(x0.^2+y0.^2+z0.^2);
% MSD_sm=smooth(t0,ind_MSD{i},0.2);
% deriv_msd_s=diff(MSD_sm)/exposure;
% D_inst=abs(deriv_msd_s/(6*10^6));
% 
% % subplot(1,3,2);hold on;
% % plot(t0,ind_MSD{i});
% subplot(1,3,2);hold on;
% plot(t0,MSD_sm);
% subplot(1,3,3);hold on;
% plot(t0(1:end-1)+0.5*exposure,D_inst);

