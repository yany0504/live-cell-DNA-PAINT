%% diffusion vs. distance (NN synapse to center of mass of the trajectory)

fileName1=[outputpath 'result_tracking_' outputname '.txt'];
fileName2=[outputpath 'Synapses1_' outputname '.mat'];
fileName3=[outputpath 'diffusion_coefficients_' outputname '.mat'];
result_tracking=textread(fileName1);
S_Synapses1=load(fileName2);
Synapses1=S_Synapses1.Synapses1;
load(fileName3);

x=result_tracking(:,1);
y=result_tracking(:,2);
z=result_tracking(:,3);
frame=result_tracking(:,4);
particle=result_tracking(:,5);

number_traj=length(Dif_track_W_ID);

% link diffusion coefficient and trajectory, then
% find the center of trajectory

for n=1:number_traj;
    dp=find(particle==Dif_track_W_ID(n,2));
    dpmax=max(dp);
    dpmin=min(dp);
    m(n)=n; %test k value
    x1=x(dpmin:dpmax);
    y1=y(dpmin:dpmax);
    z1=z(dpmin:dpmax);
    pt1=particle(dpmin:dpmax);
    %trace=[x1,y1,z1];
    frame1=frame(dpmin:dpmax);
    traces{n}=[x1,y1,z1,frame1,pt1];% creating a cell array for particle# k traces{k}
    traces_center(n,1)=mean(x1);
    traces_center(n,2)=mean(y1);
    traces_center(n,3)=mean(z1);
    %traces_center(n,4)=Dif_track_W_ID(n,2);
end

% find the distance between the center of the
% trajectories and 
nSyn=length(Synapses1);
for rg=1:nSyn
    syn_rg=Synapses1{rg}; % load data
    syn_length=length(syn_rg);
    % Find center of the synapse
    Ctr0=mean(syn_rg);  % center of synapse(rg)
    synapse_center(rg,:)=Ctr0(1,1:3);        
end
distance=zeros(number_traj,1);

for k=1:number_traj;
    for j=1:nSyn
        dist_temp=(traces_center(k,:)-synapse_center(j,:)).^2;
        dist1(j,1)=sqrt(dist_temp(1,1)+dist_temp(1,2)+dist_temp(1,3));
    end
    distance(k,1)=min(dist1);
end

distance=distance/1000; % unit conversion to um.
output_diff_dist=[Dif_track_W_ID(:,1) log10(Dif_track_W_ID(:,1)) distance Dif_track_W_ID(:,2)];
xlswrite(strcat(outputpath,['DiffusionCoeff_Distance_' [outputname '_' outputindex{1}] '.xlsx']), output_diff_dist);

