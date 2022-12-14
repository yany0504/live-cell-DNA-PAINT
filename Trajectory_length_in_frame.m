File='E:\Data\DNA-PAINT(Live-cell)\SPT_tracking\Basal_tracking\210304_DNA-PAINT_1xR3_antiHomer1-mGeos_3D_DIV15_HP\7xR3_1nM_R3_9nt-LD655\Results_timelength\MLE_timeadded_synaptic_Diff_Traj_Dist.xlsx';
data=xlsread(File);
exposure=0.1;
immobile=data(data(:,1)<-2.5,:);
mobile=data(data(:,1)>=-2.5,:);
t_mobile=exposure*mobile(:,5);
t_immobile=exposure*immobile(:,5);
mean(t_mobile)
std(t_mobile)
mean(t_immobile)
std(t_immobile)