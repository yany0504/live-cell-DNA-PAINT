% This script detects nearby traces based on the center of the current
% synapse.
% Need to load Synapses and run trace_center.m before this script.
clear syn_traces trace_n_syn1 syn_w_trace1 k kl mk mk_l Ctr0 Ctr00 Sig0 Sig00 nSyn nSyn2 radii rg syn_dis syn_length syn_rg syn_sub 

%===================================
% Load synapse
%load('Synapses.mat')
nSyn=length(Synapses);
% Load trace and find the centers of the traces
%load('traces_center.mat')
%===========================
% Define a range around synapse center to find traces that is close to the
% synapse
syn_dis=5000; % unit in nm
%===========================
% Loop
for rg=1:nSyn
syn_rg=Synapses{rg}; % load data
syn_length=length(syn_rg);
% Find center of the synapse
clear trace_nearby syn_traces
radii=1;
[Ctr00,Sig00]=subclust(syn_rg,radii); % Find center of a subcluster and the range of influence of the center
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
syn_sub(rg,:)=Ctr0(1:3);
% Find traces close to the synapse
syn_traces(:,1)=(traces_center(:,1)-Ctr0(1)).*(traces_center(:,1)-Ctr0(1));
syn_traces(:,2)=(traces_center(:,2)-Ctr0(2)).*(traces_center(:,2)-Ctr0(2));
syn_traces(:,3)=(traces_center(:,3)-Ctr0(3)).*(traces_center(:,3)-Ctr0(3));
syn_traces(:,4)=sqrt(syn_traces(:,1)+syn_traces(:,2)+syn_traces(:,3));

trace_nearby=find(syn_traces(:,4)<syn_dis);
syn_nearby{rg}=trace_nearby; 
end
syn_nearby0=find(~cellfun(@isempty,syn_nearby)); % Synapses ID with traces around 

% Plot
nSyn2=length(syn_nearby0); 
figure
hold all
% for rg=1:nSyn
%     syn_rg=Synapses{rg};
%     scatter3(syn_rg(:,1), syn_rg(:,2), syn_rg(:,3),30,'filled','c')
% end
%% Plot all tracse in the figure
clear trace_K traceID
traceID_length=size(traces')
for traceID=1:traceID_length
    trace_K=traces{traceID};
    plot3(trace_K(:,1),trace_K(:,2),trace_K(:,3),'g')
end
%clear trace_K traceID

for k=1:nSyn2
    mk=syn_nearby{syn_nearby0(k)}; % syn_nearby
    mk_l=length(mk);
    for kl=1:mk_l
        mk_trace=traces{mk(kl)};
        plot3(mk_trace(:,1),mk_trace(:,2),mk_trace(:,3),'r');
    end
    mk_syn=Synapses{syn_nearby0(k)};
    scatter3(mk_syn(:,1), mk_syn(:,2), mk_syn(:,3),10,'filled','b');
    %% Output cell array of trace that is close to synapses
    mk_trace1=mk_trace;
    mk_trace1(:,6)=syn_nearby0(k)*ones(length(mk_trace),1);
    trace_n_syn1{k}=mk_trace1;% for each cell, x,y,z,frame,particle,synapse
    %% Output cell array of synapses with traces around
    syn_w_trace1{k}=mk_syn; % for each cell, x,y,z,class,synapse
    %% Clear 
    clear mk_1 mk_trace mk_syn mk_trace1
end


%trace_n_syn=trace_n_syn1';
%syn_w_trace=syn_w_trace1';
%mTraceID=cell2mat(syn_nearby');
%TraceID=unique(mTraceID);

hold off 
axis equal
clear syn_traces trace_n_syn1 syn_w_trace1 k kl mk mk_l Ctr0 Ctr00 Sig0 Sig00 nSyn nSyn2 radii rg syn_dis syn_length syn_rg syn_sub trace_nearby



