% Syn_COAR_all: For all the synapses, syanpse center, orientation, axis, range
% This script rotates a 3D cluster data (synapse) to its principal axis and
% calculates its center and sigma. A synaptic range is decided based on the
% ellipsoid fitted to the original scattered data. Than the synapse data
% and ellipsoid is rotated back to the original coordinate system. 

% Get data
% load('Synapses.mat')
syn_n=length(Synapses);
%rg=1;
figure
hold all
for rg=1:syn_n % Define which synapse
%% load data
syn_rg=Synapses{rg}; % load data
if isempty(syn_nearby{rg})==0 % find traces near this synapse
    syn_trace_rg=syn_nearby{rg}; % syn_trace_rg is the ID of traces that is near synapse rg.
else
    syn_trace_rg=0;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Principal analysis
% Find the principal axis. Keep the scattered data centeted, so that the
% new scattered data from PCA can be rotated easily back to the original
% coordinate system. 
[SynR_coeff, SynR_score, latent]=pca(syn_rg(:,1:3),'Centered',false); 

syn_rr=SynR_score*inv(SynR_coeff); % syn_rr is the data rotated back to the original coordinate system.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Test if the transformation is correct
% figure;
% hold on
% scatter3(syn_rr(:,1),syn_rr(:,2),syn_rr(:,3))
% scatter3(syn_rg(:,1),syn_rg(:,2),syn_rg(:,3),5,'filled','r')
% hold off

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Find the center of synapse
radii=1;
[Ctr11,Sig11] = subclust(SynR_score,radii); % find center of a subcluster and the range of influence of the center in PCA coordinates
%% make sure Ctr0 has only 1 row
if length(Ctr11(:,1))==1
    Ctr1=Ctr11;
else
    Ctr1=mean(Ctr11);
end

if length(Sig11(:,1)==1)
    Sig1=Sig11;
else
    Sig1=mean(Sig11);
end
[Ctr00,Sig00] = subclust(syn_rg,radii); %   find center of a subcluster and the range of influence of the center
%% make sure Ctr0 has only 1 row
if length(Ctr00(:,1))==1
    Ctr0=Ctr00;
else
    Ctr0=mean(Ctr00);
end
if length(Sig00(:,1)==1)
    Sig0=Sig00;
else
    Sig0=mean(Sig00);
end
%% plot cluster
scatter3(syn_rg(:,1), syn_rg(:,2), syn_rg(:,3),30,'filled','b'); % Plot scattered data in original space


% plot center
scatter3(Ctr0(1), Ctr0(2),Ctr0(3), 100, 'filled','r'); % Plot center of cluster

% Pricipal axis
axis1_x1=-SynR_coeff(1,1)*500+Ctr0(1);
axis1_x2=SynR_coeff(1,1)*500+Ctr0(1);
axis1_y1=-SynR_coeff(2,1)*500+Ctr0(2);
axis1_y2=SynR_coeff(2,1)*500+Ctr0(2);
axis1_z1=-SynR_coeff(3,1)*500+Ctr0(3);
axis1_z2=SynR_coeff(3,1)*500+Ctr0(3);
plot3([axis1_x1 axis1_x2], [axis1_y1 axis1_y2], [axis1_z1 axis1_z2],'m','LineWidth',2)
% Second axis
axis2_x1=-SynR_coeff(1,2)*500+Ctr0(1);
axis2_x2=SynR_coeff(1,2)*500+Ctr0(1);
axis2_y1=-SynR_coeff(2,2)*500+Ctr0(2);
axis2_y2=SynR_coeff(2,2)*500+Ctr0(2);
axis2_z1=-SynR_coeff(3,2)*500+Ctr0(3);
axis2_z2=SynR_coeff(3,2)*500+Ctr0(3);
plot3([axis2_x1 axis2_x2], [axis2_y1 axis2_y2], [axis2_z1 axis2_z2],'c','LineWidth',2)
% Third axis, this is the main axis of a synapse, i.e. the synaptic
% orientation
axis3_x1=-SynR_coeff(1,3)*500+Ctr0(1);
axis3_x2=SynR_coeff(1,3)*500+Ctr0(1);
axis3_y1=-SynR_coeff(2,3)*500+Ctr0(2);
axis3_y2=SynR_coeff(2,3)*500+Ctr0(2);
axis3_z1=-SynR_coeff(3,3)*500+Ctr0(3);
axis3_z2=SynR_coeff(3,3)*500+Ctr0(3);
plot3([axis3_x1 axis3_x2], [axis3_y1 axis3_y2], [axis3_z1 axis3_z2],'k','LineWidth',2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot ellipsoid
[Ex,Ey,Ez] = ellipsoid(Ctr1(1),Ctr1(2),Ctr1(3),Sig1(1),Sig1(2),Sig1(3),20); 
%surf(Ex,Ey,Ez)
%Reshape the grid for surface to make a matrix in X Y and Z for the
%ellipsoid
Ex_X=reshape(Ex,21*21,1);
Ey_Y=reshape(Ey,21*21,1);
Ez_Z=reshape(Ez,21*21,1);
EXYZ=[Ex_X,Ey_Y,Ez_Z];
% Rotate the ellipsoid back to in the original space
E_back=EXYZ*inv(SynR_coeff);
Ex_X_b=reshape(E_back(:,1),21,21);
Ey_Y_b=reshape(E_back(:,2),21,21);
Ez_Z_b=reshape(E_back(:,3),21,21);
surf(Ex_X_b,Ey_Y_b,Ez_Z_b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create another ellipsoid with defined semi-axis such that within the
% ellisoid, the data is considered synaptic and outside of it, the data is
% considered extra-synaptic.
Syn_buff=1000; % The buffer layer outside synapse
[Syn_Ex, Syn_Ey, Syn_Ez]=ellipsoid(Ctr1(1),Ctr1(2),Ctr1(3),Sig1(1)+Syn_buff,Sig1(2)+Syn_buff,Sig1(3)+Syn_buff,20);

Syn_Ex_X=reshape(Syn_Ex,21*21,1);
Syn_Ey_Y=reshape(Syn_Ey,21*21,1);
Syn_Ez_Z=reshape(Syn_Ez,21*21,1);
Syn_EXYZ=[Syn_Ex_X,Syn_Ey_Y,Syn_Ez_Z];
% Rotate the ellipsoid back to in the original space
Syn_E_back=Syn_EXYZ*inv(SynR_coeff);

Syn_Ex_b=reshape(Syn_E_back(:,1),21,21);
Syn_Ey_b=reshape(Syn_E_back(:,2),21,21);
Syn_Ez_b=reshape(Syn_E_back(:,3),21,21);
surf(Syn_Ex_b, Syn_Ey_b, Syn_Ez_b);

alpha(.1)


%% In and out of synapse
clear k clear mTraces_rg traces2
% Get coordinate of traces near synapse_rg 
if isempty(syn_nearby{rg})==0
    for k=1:length(syn_trace_rg)
        traces2{k}=traces{syn_trace_rg(k)}; 
    end
    mTraces_rg1=cell2mat(traces2');
    mTraces_rg=mTraces_rg1(:,1:3);
else mTraces_rg=zeros(3);
end
      
IN=inhull(mTraces_rg,Syn_E_back);   
IN_ptc0=[IN, IN, IN].*mTraces_rg;
IN_ptc=IN_ptc0(find(IN_ptc0(:,1)),:);
OUT_ptc0=mTraces_rg-IN_ptc0;
OUT_ptc=OUT_ptc0(find(OUT_ptc0(:,1)),:);
scatter3(IN_ptc(:,1),IN_ptc(:,2),IN_ptc(:,3),30,'filled','r')
scatter3(OUT_ptc(:,1),OUT_ptc(:,2),OUT_ptc(:,3),30,'filled','g')


% Output
clear IN_trace0 IN_trace1 
IN_trace0=[IN,IN,IN,IN,IN].*mTraces_rg1;
IN_trace1=IN_trace0(find(IN_trace0(:,1)),:);
if isempty(IN_trace1)==0
    IN_trace{rg}=IN_trace1;
else
    IN_trace{rg}=zeros(1,5);
end

OUT_trace0=mTraces_rg1-IN_trace0;
OUT_trace1=OUT_trace0(find(OUT_trace0(:,1)),:);
if isempty(OUT_trace1)==0
    OUT_trace{rg}=OUT_trace1;
else
    OUT_trace{rg}=zeros(1,5);
end
% % Test
% 
% ptc=-500+1000*rand(100,3);
% ptc_ctr=[ptc(:,1)+Ctr0(1), ptc(:,2)+Ctr0(2), ptc(:,3)+Ctr0(3)]
% IN=inhull(ptc_ctr,Syn_E_back);
% IN_ptc0=[IN, IN, IN].*ptc_ctr;
% IN_ptc=IN_ptc0(find(IN_ptc0(:,1)),:);
% OUT_ptc0=ptc_ctr-IN_ptc0;
% OUT_ptc=OUT_ptc0(find(OUT_ptc0(:,1)),:);
% 
% scatter3(IN_ptc(:,1),IN_ptc(:,2),IN_ptc(:,3),30,'filled','r')
% scatter3(OUT_ptc(:,1),OUT_ptc(:,2),OUT_ptc(:,3),30,'filled','g')
%% %%
grid off
axis equal

end

hold off

clear Ctr0 Ctr00 Ctr1 Ctr11 EXYX E_back Ex Ex_X Ex_X_b Ey Ey_Y Ey_Y_b Ez Ez_Z Ez_Z_b 
clear Sig0 Sig00 Sig1 Sig11 SynR_coeff SynR_score Syn_EXYZ Syn_E_back Syn_Ex Syn_Ex_X Syn_Ex_b Syn_Ey Syn_Ey_Y Syn_Ey_b
clear Syn_Ez Syn_Ez_Z Syn_Ez_b Syn_buff k latentmTraces_rg mTraces_rg1 radii rg syn_n syn_rg syn_rr
clear axis1_x1 axis1_x2 axis1_y1 axis1_y2 axis1_z1 axis1_z2 
clear axis2_x1 axis2_x2 axis2_y1 axis2_y2 axis2_z1 axis2_z2
clear axis3_x1 axis3_x2 axis3_y1 axis3_y2 axis3_z1 axis3_z2