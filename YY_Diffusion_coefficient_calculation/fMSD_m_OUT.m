% This script calculates the MSD of traces with give array of specific
% particle ID.
%% Get trace ID
% Must run nearby_trace.m before running this script.
% clear mTraceIN0 mTraceIN TraceIN
% mTraceIN0=cell2mat(IN_trace');
% mTraceIN=mTraceIN0(find(mTraceIN0(:,5)),:);
% TraceIN=unique(mTraceIN(:,5)); % get the unique particle ID for the nearby traces and arrange in acsend order.
% %%

%%
%This m_file seperates particles according to their ID in the 'result', and
%calculates diffusion coeficient for each particle by fitting straight
%lines to the first 4 points of MSD/dt of the particle.
%of the tracking file and outputs the submatrices, e.g. sub_k is the
%submatrix of particle# k.
tSteps=10;
%% load file
%make variable x y z frame number and particle ID
x=mTraceOUT(:,1);
y=mTraceOUT(:,2);
z=mTraceOUT(:,3); % 0.79 is the correction factor for RI mismatch for water/oil for an 1.40 NA objective
frame=mTraceOUT(:,4);
particle=mTraceOUT(:,5);
d=length(particle); % total row number
TraceOUT=unique(mTraceOUT(:,5));
ID_length=length(TraceOUT); % length of particle number
%% MSD calculation
%% MSD calculation
for k1=1:ID_length
        k=TraceOUT(k1);
        dp=find(particle==k); % dimension of matrix of particle# K 
        dplength=length(dp);
        if dplength<2
            Dif(k,:)=0;
        else
        dpmax=max(dp); % for particleK, the first row number
        dpmin=min(dp); % for particleK, the last row number
        m(k)=k % test k value
        x1=x(dpmin:dpmax);
        y1=y(dpmin:dpmax);
        z1=z(dpmin:dpmax);
        frame1=frame(dpmin:dpmax);
        particle1=particle(dpmin:dpmax);
        Trace_OUT{k}=[x1,y1,z1,frame1,particle1];
        % eval(['sub_' num2str(k) '=[x1,y1,z1,frame1]']);% creating submatrix for particle# k sub_k 
        [MSD00,d2r0,counts]=fMSD(x1,y1,z1,frame1,dpmax,dpmin,tSteps); % Call function MSD, MSD.m file must be in the same folder as this file
        cutoff=10;  %cutoff for data points 
        ind=find(counts>cutoff-1); % find index of counts that is equate to and above cutoff
        MSDCF=MSD00(ind); % find the MSD for those index from last line
        %ind1=[0,ind]; % Add (0,0) as the first point of the curve
        %MSDCF1=[0,MSDCF]; % Add (0,0) as the first point of the curve
        % ind1=ind(1:4);
        % MSDCF1=MSDCF(1:4);
        indlength=length(ind);
        if indlength>=4
            ind1=ind(1:4);
            MSDCF1=MSDCF(1:4);
            f = fit(ind1',MSDCF1','poly1');
            Dif(k,:)=coeffvalues(f); % First column is the diffusion coeff (slope). To get um^2/s multiply by 2*10^(-5), to get cm^2/s, multiply by 2*10^(-13)
        else Dif(k,:)=0;
        end
        end
        plot(ind1,MSDCF1)
        hold all
        %eval(['MSD00_' num2str(k) '=[MSD00]']);
        %eval(['d2r0_' num2str(k) '=[d2r0]']);
        clear dp dplength dpmax dpmin x1 y1 z1 frame1
end 
hold off
dt=0.05; %Default time step 0.05s=50ms
Dif1=Dif(:,1);
Dif_positive=Dif1(find(Dif1>0));
Dif_time_OUT=Dif_positive/(dt*2*3); %D for time step dt, for 3D
Dif_track_OUT=Dif_positive/(dt*2*3*10^6);
hist(Dif_track_OUT, 15);
clear x y z frame particle k d m x1 y1 z1 frame1 dp Dif1 Dif_positive


