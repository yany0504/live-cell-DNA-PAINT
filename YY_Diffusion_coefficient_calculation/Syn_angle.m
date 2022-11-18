%   This script calculates the angle alpha between z axis and the major
%   axis of synapses

syn_n=length(Synapses);
for rg=1:syn_n
syn_rg=Synapses{rg}; % load data
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Principal analysis
% Find the principal axis. Keep the scattered data centeted, so that the
% new scattered data from PCA can be rotated easily back to the original
% coordinate system. 
[SynR_coeff, SynR_score, latent]=pca(syn_rg(:,1:3),'Centered',false); 

syn_rr=SynR_score*inv(SynR_coeff); % syn_rr is the data rotated back to the original coordinate system.
%%  Calculate the angle between z-axis and the major axis of the synapse
Axis_z=[0 0 1];
norm_coeff=sqrt(SynR_coeff(1,3)^2+SynR_coeff(2,3)^2+SynR_coeff(3,3)^2); % normalize coeff
cos_a=Axis_z*SynR_coeff(:,3)/norm_coeff;
syn_alpha(rg)=acosd(abs(cos_a)); % arccos in degrees
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%