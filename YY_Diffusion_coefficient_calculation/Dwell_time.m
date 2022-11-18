function dwell=Dwell_time(result_tracking,nearby_trace,Ctr1,Sig1,Syn_buff,rg)
%%%%%%%%%%%%%
%%%%%%trace_center
% This script finds the center of traces abtained from tracking.m
%===============================
% load data from tracking.m

x=result_tracking(:,1);
y=result_tracking(:,2);
z=0.79*result_tracking(:,3);
frame=result_tracking(:,4);
particle=result_tracking(:,5);
nearby_trace0=find(~isempty(nearby_trace));
nTrace=length(nearby_trace0);
dwell{rg}=[];
for h=1:nTrace
    dp=find(particle==nearby_trace0(h)); % dimension of matrix of particle# K
    dpmax=max(dp); % for particleK, the first row number
    dpmin=min(dp); % for particleK, the last row number
    x1=x(dpmin:dpmax);
    y1=y(dpmin:dpmax);
    z1=z(dpmin:dpmax);
    frame1=frame(dpmin:dpmax);
    particle1=particle(dpmin:dpmax);
    inSyn=find((x1-Ctr1(1)).^2/(Sig1(1)+Syn_buff)^2+(y1-Ctr1(2)).^2/(Sig1(2)+Syn_buff)^2+(z1-Ctr1(3)).^2/(Sig1(3)+Syn_buff)^2<1);
    head=1;
    tail=1;
    while head<length(inSyn)& tail<length(inSyn)
        tail=tail+1;
        if (inSyn(tail)-inSyn(tail-1))~=1
            if tail-head>8
                dwell{rg}=vertcat(dwell{rg},[frame1(tail-1)-frame1(head) particle1(head)]);
                head=tail;
            end
        end
    end
    if tail-head>8
        dwell{rg}=vertcat(dwell{rg},[frame1(tail-1)-frame1(head) particle1(head)]);
    end
end