%% Fit a 2D gaussian function to data
%% PURPOSE:  Fit a 2D gaussian centroid to simulated data
% Uses lsqcurvefit to fit
%
% INPUT:
% 
%   MdataSize: Size of nxn data matrix
%   x0 = [Amp,x0,wx,y0,wy,fi]: Inital guess parameters
%   x = [Amp,x0,wx,y0,wy,fi]: simulated centroid parameters
%   noise: noise in % of centroid peak value (x(1)
%   InterpolationMethod = 'nearest' or 'linear' or 'spline' or 'cubic'
%       used to calculate cross-section along minor and major axis
%     
%
%
% NOTE:
%   The initial values in x0 must be close to x in order for the fit
%   to converge to the values of x (especially if noise is added)
%
% OUTPUT:  non
%
% CREATED: G. Nootz  May 2012
% 
%  Modifications:
%  non
%% ---------User Input---------------------
clear
clc;
MdataSize = 50; % Size of nxn data matrix
% parameters are: [Amplitude, x0, sigmax, y0, sigmay, angel(in rad)]
x0 = [1,-2,5,3,5,0,2,3,5,4,5,0]; %Inital guess parameters
x = [2,-2.2,3,-3.4,2,+0.02*2*pi,2,3.2,3,5.4,2,+0]; %centroid parameters
noise = 10; % noise in % of centroid peak value (x(1))
InterpolationMethod = 'nearest'; % 'nearest','linear','spline','cubic'
FitForOrientation = 0; % 0: fit for orientation. 1: do not fit for orientation
%% ---Generate centroid to be fitted--------------------------------------
xin = x; 
noise = noise/100 * x(1);
[X,Y] = meshgrid(-MdataSize/2:MdataSize/2);
xdata = zeros(size(X,1),size(Y,2),2);
xdata(:,:,1) = X;
xdata(:,:,2) = Y;
[Xhr,Yhr] = meshgrid(linspace(-MdataSize/2,MdataSize/2,300)); % generate high res grid for plot
xdatahr = zeros(300,300,2);
xdatahr(:,:,1) = Xhr;
xdatahr(:,:,2) = Yhr;
%---Generate noisy centroid---------------------
Z = Two2DGaussFunctionRot(x,xdata);
Z = Z + noise*(rand(size(X,1),size(Y,2))-0.5);
%% --- Fit---------------------
if FitForOrientation == 0
    % define lower and upper bounds [Amp,xo,wx,yo,wy,fi]
    lb = [0,-MdataSize/2,0,-MdataSize/2,0,-pi/4,0,-MdataSize/2,0,-MdataSize/2,0,-pi/4];
    ub = [realmax('double'),MdataSize/2,(MdataSize/2)^2,MdataSize/2,(MdataSize/2)^2,pi/4,realmax('double'),MdataSize/2,(MdataSize/2)^2,MdataSize/2,(MdataSize/2)^2,pi/4];
    [x,resnorm,residual,exitflag] = lsqcurvefit(@Two2DGaussFunctionRot,x0,xdata,Z,lb,ub);
else
    x0 =x0(1:5);
    xin(6) = 0; 
    x =x(1:5);
    lb = [0,-MdataSize/2,0,-MdataSize/2,0];
    ub = [realmax('double'),MdataSize/2,(MdataSize/2)^2,MdataSize/2,(MdataSize/2)^2];
    [x,resnorm,residual,exitflag] = lsqcurvefit(@Two2DGaussFunctionRot,x0,xdata,Z,lb,ub);
    x(6) = 0;
end
%% ---------Plot 3D Image-------------
figure(1)
C = del2(Z);
mesh(X,Y,Z,C) %plot data
hold on
surface(Xhr,Yhr,Two2DGaussFunctionRot(x,xdatahr),'EdgeColor','none') %plot fit
axis([-MdataSize/2-0.5 MdataSize/2+0.5 -MdataSize/2-0.5 MdataSize/2+0.5 -noise noise+x(1)])
alpha(0.2)  
hold off
%% -----Plot profiles----------------
hf2 = figure(2);
set(hf2, 'Position', [20 20 950 900])
alpha(0)
subplot(4,4, [5,6,7,9,10,11,13,14,15])
imagesc(X(1,:),Y(:,1)',Z)
set(gca,'YDir','reverse')
colormap('jet')
string1 = ['       Amplitude1','    X-Coordinate1', '    X-Width1','    Y-Coordinate1','    Y-Width1','     Angle1','       Amplitude2','    X-Coordinate2', '    X-Width2','    Y-Coordinate2','    Y-Width2','     Angle2'];
string2 = ['Set     ',num2str(xin(1), '% 100.3f'),'             ',num2str(xin(2), '% 100.3f'),'         ',num2str(xin(3), '% 100.3f'),'         ',num2str(xin(4), '% 100.3f'),'        ',num2str(xin(5), '% 100.3f'),'     ',num2str(xin(6), '% 100.3f'),'     ',num2str(xin(7), '% 100.3f'),'             ',num2str(xin(8), '% 100.3f'),'         ',num2str(xin(9), '% 100.3f'),'         ',num2str(xin(10), '% 100.3f'),'        ',num2str(xin(11), '% 100.3f'),'     ',num2str(xin(12), '% 100.3f')];
string3 = ['Fit      ',num2str(x(1), '% 100.3f'),'             ',num2str(x(2), '% 100.3f'),'         ',num2str(x(3), '% 100.3f'),'         ',num2str(x(4), '% 100.3f'),'        ',num2str(x(5), '% 100.3f'),'     ',num2str(x(6), '% 100.3f'),'     ',num2str(x(7), '% 100.3f'),'             ',num2str(x(8), '% 100.3f'),'         ',num2str(x(9), '% 100.3f'),'         ',num2str(x(10), '% 100.3f'),'        ',num2str(x(11), '% 100.3f'),'     ',num2str(x(12), '% 100.3f')];
text(-MdataSize/2*0.9,+MdataSize/2*1.15,string1,'Color','red')
text(-MdataSize/2*0.9,+MdataSize/2*1.2,string2,'Color','red')
text(-MdataSize/2*0.9,+MdataSize/2*1.25,string3,'Color','red')
%% -----Calculate cross sections-------------
% generate points along horizontal axis
m = -tan(x(6));% Point slope formula
b = (-m*x(2) + x(4));
xvh = -MdataSize/2:MdataSize/2;
yvh = xvh*m + b;
hPoints = interp2(X,Y,Z,xvh,yvh,InterpolationMethod);
% generate points along vertical axis
mrot = -m;
brot = (mrot*x(4) - x(2));
yvv = -MdataSize/2:MdataSize/2;
xvv = yvv*mrot - brot;
vPoints = interp2(X,Y,Z,xvv,yvv,InterpolationMethod);
hold on % Indicate major and minor axis on plot
% % plot pints 
% plot(xvh,yvh,'r.') 
% plot(xvv,yvv,'g.')
% plot lins 
plot([xvh(1) xvh(size(xvh))],[yvh(1) yvh(size(yvh))],'r') 
plot([xvv(1) xvv(size(xvv))],[yvv(1) yvv(size(yvv))],'g') 

m2 = tan(x(12));
b2=(-m2*x(8)+x(10));
xvh2= -MdataSize/2:MdataSize/2;
yvh2 = xvh2*m2 + b2;
hPoints2 = interp2(X,Y,Z,xvh2,yvh2,InterpolationMethod);

mrot2 = -m2;
brot2 = (mrot2*x(10) - x(8));
yvv2 = -MdataSize/2:MdataSize/2;
xvv2 = yvv2*mrot2 - brot2;
vPoints2 = interp2(X,Y,Z,xvv2,yvv2,InterpolationMethod);
plot([xvh2(1) xvh2(size(xvh2))],[yvh2(1) yvh2(size(yvh2))],'r') 
plot([xvv2(1) xvv2(size(xvv2))],[yvv2(1) yvv2(size(yvv2))],'g') 
hold off
axis([-MdataSize/2-0.5 MdataSize/2+0.5 -MdataSize/2-0.5 MdataSize/2+0.5])
hPoints_tot=hPoints+hPoints2;
vPoints_tot=vPoints+vPoints2;
%%
ymin = - noise * x(1);
ymax = x(1)*(1+noise);
xdatafit = linspace(-MdataSize/2-0.5,MdataSize/2+0.5,300);
hdatafit = x(1)*exp(-(xdatafit-x(2)).^2/(2*x(3)^2))+x(7)*exp(-(xdatafit-x(8)).^2/(2*x(9)^2));
vdatafit = x(1)*exp(-(xdatafit-x(4)).^2/(2*x(5)^2))+x(7)*exp(-(xdatafit-x(10)).^2/(2*x(11)^2));
subplot(4,4, [1:3])
xposh = (xvh-x(2))/cos(x(6))+x(2);% correct for the longer diagonal if fi~=0
plot(xposh,hPoints_tot,'r.',xdatafit,hdatafit,'black')
axis([-MdataSize/2-0.5 MdataSize/2+0.5 ymin*1.1 ymax*1.1])
subplot(4,4,[8,12,16])
xposv = (yvv-x(4))/cos(x(6))+x(4);% correct for the longer diagonal if fi~=0
plot(vPoints_tot,xposv,'g.',vdatafit,xdatafit,'black')
axis([ymin*1.1 ymax*1.1 -MdataSize/2-0.5 MdataSize/2+0.5])
set(gca,'YDir','reverse')
figure(gcf) % bring current figure to front