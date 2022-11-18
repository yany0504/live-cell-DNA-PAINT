% This script output matrix of traces with give array of specific
% particle ID for all traces near synapses and amoung these traces, the
% traces (or part of traces that is) in the synapse and out of synapse. 
% This script need input from nearby_traces.m syn_COAR.m
%% Get trace ID
% Must run nearby_trace.m before running this script.
mTraceID=cell2mat(syn_nearby');
TraceID=unique(mTraceID); % get the unique particle ID for the nearby traces and arrange in acsend order.



%% load file
%make variable x y z frame number and particle ID
x=result(:,1);
y=result(:,2);
z=0.79*result(:,3); % 0.79 is the correction factor for RI mismatch for water/oil for an 1.40 NA objective
frame=result(:,4);
particle=result(:,5);
d=length(particle); % total row number
ID_length=length(TraceID); % length of particle number
%% MSD calculation
for k1=1:ID_length
        k=TraceID(k1);
        dp=find(particle==k); % dimension of matrix of particle# K 
        dplength=length(dp);
        dpmax=max(dp); % for particleK, the first row number
        dpmin=min(dp); % for particleK, the last row number
        m(k)=k; % test k value
        x1=x(dpmin:dpmax);
        y1=y(dpmin:dpmax);
        z1=z(dpmin:dpmax);
        frame1=frame(dpmin:dpmax);
        particle1=particle(dpmin:dpmax);
        Trace{k}=[x1,y1,z1,frame1,particle1];
        clear dp dplength dpmax dpmin x1 y1 z1 frame1
end 
mTrace=cell2mat(Trace');
mTraceIN0=cell2mat(IN_trace');
mTraceIN1=mTraceIN0(find(mTraceIN0(:,5)),:);
mTraceIN=intersect(mTraceIN1, mTraceIN1,'rows','stable');
mTraceOUT=setdiff(mTrace,mTraceIN,'rows','stable');

clear x y z frame particle k d m x1 y1 z1 frame1 particle1 dp ID_length 


