function [ MSD0, d2r, counts ] = fMSD_vect( x1,y1,z1,frame1,dpmax,dpmin,tSteps )
%fMSD function calculates the MSD of a 3D diffusing particle
%   Input: x1,y1,z1,frame1 are from the result of particle tracking matrix,
%   this function has to work with MSD_m because of the input. tSteps is
%   the biggest number of steps in time for MSD. E.g: tSteps=10: Calculate MSD for the
%   first 10 time steps instead of calculating MSD for all the time steps.
%   dpmax: for particleK, the last row number
%   dpmin: for particleK, the first row number
%   Output: MSD0 is a column vector MSD0(c) for mean square diviation at
%   timestep c. d2r is a matrix for distance square. d2r (a,c) is the ath
%   value of distance square at timestep c. Counts is a column vector for
%   number of points for timestep c.

x0=x1;
y0=y1;
z0=z1;
frame0=frame1;
dpmax0=dpmax;
dpmin0=dpmin;
MSD0=zeros(1,max(frame0)-min(frame0));
if tSteps<=(max(frame0)-min(frame0))
    d2r=zeros(dpmax-dpmin,tSteps);
    for c=1:tSteps % c: time steps.  from 1 to tSteps % previous c=1:max(frame0)-min(frame0)
        
        [B, A]=meshgrid(2:(dpmax0-dpmin0+1),1:(dpmax0-dpmin0));
        [a,b]=find((frame0(B)-frame0(A))==c);
        b=b+1;
        %                      frameab(a,b)=frame0(b)-frame0(a);
        d2r(a,c)=(x0(b)-x0(a)).^2+(y0(b)-y0(a)).^2+(z0(b)-z0(a)).^2; % SD for frame difference c
        % For different b, d2r(a,c) remains the same. Because: for fixed c and fixed a, frame0(b)=c+frame0(a) is fixed. So b is fixed.
        %eval(['deltaRsq_' num2str(c) '=d2r'])%deltaR^2 for each value c
        MSD0(c)=mean(nonzeros(d2r(:,c)));
        counts(c)=length(nonzeros(d2r(:,c))); % counts for each timestep C
        %cutoff=10;  %cutoff for data points
        %ind=find(counts>cutoff-1) % find index of counts that is equate to and above cutoff
        %MSDCF=MSD0(ind) % find the MSD for those index from last line
    end
    
else
    
    d2r=zeros(dpmax-dpmin,max(frame0)-min(frame0));
    for c=1:max(frame0)-min(frame0) % c from 1 to biggest frame difference
        [B, A]=meshgrid(2:(dpmax0-dpmin0+1),1:(dpmax0-dpmin0));
        [a,b]=find((frame0(B)-frame0(A))==c); 
        b=b+1;
        d2r(a,c)=(x0(b)-x0(a)).^2+(y0(b)-y0(a)).^2+(z0(b)-z0(a)).^2; % SD for frame difference c
        %eval(['deltaRsq_' num2str(c) '=d2r'])%deltaR^2 for each value c        
        MSD0(c)=mean(nonzeros(d2r(:,c)));
        counts(c)=length(nonzeros(d2r(:,c))); % counts for each timestep C
        %cutoff=10;  %cutoff for data points
        %ind=find(counts>cutoff-1) % find index of counts that is equate to and above cutoff
        %MSDCF=MSD0(ind) % find the MSD for those index from last line
        
    end
end



%clear a b frameab c x0 y0 z0 frame0
