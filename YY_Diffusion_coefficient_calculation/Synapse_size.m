% This script calculates the size of synapses found by synapse_2.m. X, Y ,Z
% and volumn is calculated
tic;
[synNumb,synNull]=size(Synapses); % get number of synapses
SynapseSize=zeros(synNumb,5);
Numb=[1:synNumb];
sizeX=zeros(synNumb,1);
sizeY=zeros(synNumb,1);
sizeZ=zeros(synNumb,1);
Vol=zeros(synNumb,1);
figure
%Volume and size
for g=1:synNumb
    m=Synapses{g};
    SynX=m(:,1);
    SynY=m(:,2);
    SynZ=m(:,3);
    %Volume
    TRI = DelaunayTri(SynX,SynY,SynZ);
    [ch v] = convexHull(TRI);
    Vol(g)=v;
    K{g}=ch;
    hold all
    trisurf(ch, TRI.X(:,1),TRI.X(:,2),TRI.X(:,3))
    %Area
%     DT_XY = DelaunayTri(SynX,SynY);
%     [cha1 a1] = convexHull(DT_XY);
%     Area_XY(g)=a1;
%     DT_YZ = DelaunayTri(SynY,SynZ);
%     [cha2 a2] = convexHull(DT_YZ);
%     Area_YZ(g)=a2;
%     DT_XZ = DelaunayTri(SynX,SynZ);
%     [cha3 a3] = convexHull(DT_XZ);
%     Area_XZ(g)=a3;
    %linear size
    sizeX(g)=max(SynX)-min(SynX);
    sizeY(g)=max(SynY)-min(SynY);
    sizeZ(g)=max(SynZ)-min(SynZ);
    
    
end
SynapseSize=[Numb, Vol,sizeX, sizeY, sizeZ];
hold off
axis equal
takeTake=toc;