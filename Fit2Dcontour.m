function [x,resnorm,residual,exitflag]=Fit2Dcontour(N,c,x0,rotation);


%% ---------User Input---------------------
% parameters are: [Amplitude1, x10, sigmax1, y10, sigmay1, angle1(in rad),Amplitude2, x20, sigmax2, y20, sigmay2, angle2(in rad)]
x0 = [2,-1.5,1,-0.5,1,0.5,   1,-3,0.5,-1,0.5,0.2]; %Inital guess parameters
InterpolationMethod = 'nearest'; % 'nearest','linear','spline','cubic'
FitForOrientation = rotation; % 0: fit for orientation. 1: do not fit for orientation
%% ---Generate centroid to be fitted--------------------------------------

[X,Y] = meshgrid(c{1},c{2});
xdata = zeros(size(X,1),size(Y,2),2);
xdata(:,:,1) = X;
xdata(:,:,2) = Y;
[Xhr,Yhr] = meshgrid(linspace(c{1}(1),c{1}(end),300),linspace(c{2}(1),c{2}(end),300)); % generate high res grid for plot
xdatahr = zeros(300,300,2);
xdatahr(:,:,1) = Xhr;
xdatahr(:,:,2) = Yhr;
%---Generate noisy centroid---------------------
Z = N';
%% --- Fit---------------------
% define lower and upper bounds [Amp1,x1o,wx1,y1o,wy1,fi1,Amp2,x2o,wx2,y2o,wy2,fi2]
if FitForOrientation == 0;
    lb = [0,-5,0,-2,0,-pi/4,0,-5,0,-2,0,-pi/4];
    ub = [realmax('double'),0,5,1,5,pi/2,realmax('double'),0,5,1,5,pi/2];
    opts=optimoptions('lsqcurvefit');
    opts.MaxFunctionEvaluations = 100000;
    opts.MaxIterations = 100000;
    [x,resnorm,residual,exitflag] = lsqcurvefit(@Two2DGaussFunctionRot,x0,xdata,Z,lb,ub,opts);
else
   x0norot=[x0(1:5),x0(7:11)];
   lb = [0,-5,0,-2,0,0,-5,0,-2,0];
   ub = [realmax('double'),0,5,1,5,realmax('double'),0,5,1,5];
   opts=optimoptions('lsqcurvefit');
   opts.MaxFunctionEvaluations = 100000;
   opts.MaxIterations = 100000;
   [x,resnorm,residual,exitflag] = lsqcurvefit(@Two2DGaussFunction,x0norot,xdata,Z,lb,ub,opts);
    
end
    

%% ---------Plot 3D Image-------------
figure;
C = del2(Z);
mesh(X,Y,Z) %plot data
hold on
if FitForOrientation == 0;
    surface(Xhr,Yhr,Two2DGaussFunctionRot(x,xdatahr),'EdgeColor','none') %plot fit
else
    surface(Xhr,Yhr,Two2DGaussFunction(x,xdatahr),'EdgeColor','none') %plot fit
end

alpha(0.5)
ax=gca;
ax.ZAxis.FontSize=15;
ax.YAxis.FontSize=15;
ax.XAxis.FontSize=15;
set(gca,'FontWeight','bold');
xlabel('Diffusion coefficient (\mum^2/s, log)','fontweight','bold','FontSize',18,'Rotation',-20);
ylabel('Trajectory range (\mum, log)','fontweight','bold','FontSize',18,'Rotation',20);
zlabel('Count','fontweight','bold','FontSize',18);
view([45 20]);
set(gcf, 'Position', [100 100 1000 800]);
hold off

%% -----Plot profiles----------------
% hf2 = figure(2);
% set(hf2, 'Position', [20 20 950 900])
% alpha(0)
% subplot(4,4, [5,6,7,9,10,11,13,14,15])
% imagesc(X(1,:),Y(:,1)',Z)
% set(gca,'YDir','reverse')
% colormap('jet')
% string1 = ['       Amplitude1','    X-Coordinate1', '    X-Width1','    Y-Coordinate1','    Y-Width1','     Angle1','       Amplitude2','    X-Coordinate2', '    X-Width2','    Y-Coordinate2','    Y-Width2','     Angle2'];
% string3 = ['Fit      ',num2str(x(1), '% 100.3f'),'             ',num2str(x(2), '% 100.3f'),'         ',num2str(x(3), '% 100.3f'),'         ',num2str(x(4), '% 100.3f'),'        ',num2str(x(5), '% 100.3f'),'     ',num2str(x(6), '% 100.3f'),'     ',num2str(x(7), '% 100.3f'),'             ',num2str(x(8), '% 100.3f'),'         ',num2str(x(9), '% 100.3f'),'         ',num2str(x(10), '% 100.3f'),'        ',num2str(x(11), '% 100.3f'),'     ',num2str(x(12), '% 100.3f')];
% text(1,1,string1,'Color','red')
% text(2,2,string3,'Color','red')
%% -----Calculate cross sections-------------
% % generate points along horizontal axis
% m = -tan(x(6));% Point slope formula
% b = (-m*x(2) + x(4));
% xvh = -5:0;
% yvh = xvh*m + b;
% hPoints = interp2(X,Y,Z,xvh,yvh,InterpolationMethod);
% % generate points along vertical axis
% mrot = -m;
% brot = (mrot*x(4) - x(2));
% yvv = -1.5:0.5;
% xvv = yvv*mrot - brot;
% vPoints = interp2(X,Y,Z,xvv,yvv,InterpolationMethod);
% hold on % Indicate major and minor axis on plot
% % % plot pints 
% % plot(xvh,yvh,'r.') 
% % plot(xvv,yvv,'g.')
% % plot lins 
% plot([xvh(1) xvh(size(xvh))],[yvh(1) yvh(size(yvh))],'r') 
% plot([xvv(1) xvv(size(xvv))],[yvv(1) yvv(size(yvv))],'g') 
% 
% m2 = tan(x(12));
% b2=(-m2*x(8)+x(10));
% 
% xvh2= -5:0;
% yvh2 = xvh2*m2 + b2;
% hPoints2 = interp2(X,Y,Z,xvh2,yvh2,InterpolationMethod);
% 
% mrot2 = -m2;
% brot2 = (mrot2*x(10) - x(8));
% yvv2 = -1.5:0.5;
% xvv2 = yvv2*mrot2 - brot2;
% vPoints2 = interp2(X,Y,Z,xvv2,yvv2,InterpolationMethod);
% plot([xvh2(1) xvh2(size(xvh2))],[yvh2(1) yvh2(size(yvh2))],'r') 
% plot([xvv2(1) xvv2(size(xvv2))],[yvv2(1) yvv2(size(yvv2))],'g') 
% hold off
% %axis([-MdataSize/2-0.5 MdataSize/2+0.5 -MdataSize/2-0.5 MdataSize/2+0.5])
% hPoints_tot=hPoints+hPoints2;
% vPoints_tot=vPoints+vPoints2;
% %%
% ymin = - 1.5;
% ymax = 0.5;
% xdatafit = linspace(-5,0,300);
% xdatafit2 = linspace(-1.5,0.5,300);
% hdatafit = x(1)*exp(-(xdatafit-x(2)).^2/(2*x(3)^2))+x(7)*exp(-(xdatafit-x(8)).^2/(2*x(9)^2));
% vdatafit = x(1)*exp(-(xdatafit2-x(4)).^2/(2*x(5)^2))+x(7)*exp(-(xdatafit2-x(10)).^2/(2*x(11)^2));
% subplot(4,4, [1:3])
% xposh = (xvh-x(2))/cos(x(6))+x(2);% correct for the longer diagonal if fi~=0
% plot(xposh,hPoints_tot,'r.',xdatafit,hdatafit,'black')
% %axis([-MdataSize/2-0.5 MdataSize/2+0.5 ymin*1.1 ymax*1.1])
% subplot(4,4,[8,12,16])
% xposv = (yvv-x(4))/cos(x(6))+x(4);% correct for the longer diagonal if fi~=0
% plot(xposv,vPoints_tot,'g.',vdatafit,xdatafit2,'black')
% %axis([ymin*1.1 ymax*1.1 -MdataSize/2-0.5 MdataSize/2+0.5])
% set(gca,'YDir','reverse')
% figure(gcf) % bring current figure to front