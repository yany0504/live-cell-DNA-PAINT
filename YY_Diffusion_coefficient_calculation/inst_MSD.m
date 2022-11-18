function inst_msd= inst_MSD( x1,y1,z1,frame1)
x0=x1;
y0=y1;
z0=z1;
dt=0.05; %50ms
frame0=frame1;
span=6;%For the inst_msd at position 6, calculate from 0 to 12
head=min(find(frame0-span>=0));
tail=max(find(frame0+span<=max(frame0)));
inst_msd=zeros(1,head-tail+1);
Dif=zeros(head-tail+1,2);
n=3;
for j=head:tail
    dpmin0=min(find(frame0-frame0(j)+span>=0));
    dpmax0=max(find(frame0-frame0(j)-span<=0));
    tSteps=3;
    MSD0=zeros(1,dpmax0-dpmin0);
    d2r=zeros(dpmax0-dpmin0,tSteps);
    for c=1:tSteps % c: time steps.  from 1 to tSteps % previous c=1:max(frame0)-min(frame0)
        [B, A]=meshgrid(2:(dpmax0-dpmin0+1),1:(dpmax0-dpmin0));
        [a,b]=find((frame0(B)-frame0(A))==c);
        b=b+1;
        d2r(a,c)=(x0(b)-x0(a)).^2+(y0(b)-y0(a)).^2+(z0(b)-z0(a)).^2; % SD for frame difference c
        % For different b, d2r(a,c) remains the same. Because: for fixed c and fixed a, frame0(b)=c+frame0(a) is fixed. So b is fixed.
        %eval(['deltaRsq_' num2str(c) '=d2r'])%deltaR^2 for each value c
        MSD0(c)=mean(nonzeros(d2r(:,c)));
    end
    MSDCF1=MSD0(1:3);
    Dif(j,:)= polyfit((1:3)',MSDCF1',1);
    inst_msd(j)=Dif(j,1)/2/n/dt;
    % ORIGINAL        f = fit(ind1',MSDCF1','poly1');
    % ORIGINAL       Dif(j,:)=coeffvalues(f); % First column is the diffusion coeff (slope). To get um^2/s multiply by 2*10^(-5), to get cm^2/s, multiply by 2*10^(-13)
end
index=(1:size(inst_msd,2))';
index=index(~isnan(inst_msd));
inst_msd=inst_msd(~isnan(inst_msd));
inst_msd=[index inst_msd'];
end

%clear a b frameab c x0 y0 z0 frame0
