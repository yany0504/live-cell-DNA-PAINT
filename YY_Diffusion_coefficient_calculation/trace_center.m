% This script finds the center of traces abtained from tracking.m
%===============================
% load data from tracking.m
x=result(:,1);
y=result(:,2);
z=0.79*result(:,3);
frame=result(:,4);
particle=result(:,5);
d=length(particle); %total row number
ptotal=particle(d); %total particle number
for k=1:ptotal
        dp=find(particle==k); %dimension of matrix of particle# K 
        dplength=length(dp);
        dpmax=max(dp);
        dpmin=min(dp);
        m(k)=k; %test k value
        x1=x(dpmin:dpmax);
        y1=y(dpmin:dpmax);
        z1=z(dpmin:dpmax);
        pt1=particle(dpmin:dpmax);
        %trace=[x1,y1,z1];
        frame1=frame(dpmin:dpmax);
        traces{k}=[x1,y1,z1,frame1,pt1];% creating a cell array for particle# k traces{k}
        traces_center(k,1)=mean(x1);
        traces_center(k,2)=mean(y1);
        traces_center(k,3)=mean(z1);
end

clear x y z frame particle k d m x1 y1 z1 frame1 dp ptotal dpmax dplength dpmin pt1
