function F = Two2DGaussFunction(x,xdata)
%% x = [Amp1, x01, wx1, y01, wy1, Amp2, x02, wx2, y02, wy2]
%[X,Y] = meshgrid(x,y) 
%  xdata(:,:,1) = X
%  xdata(:,:,2) = Y           
% Mrot = [cos(fi) -sin(fi); sin(fi) cos(fi)]
%%
x1datarot(:,:,1)= xdata(:,:,1)*cos(0) - xdata(:,:,2)*sin(0);
x1datarot(:,:,2)= xdata(:,:,1)*sin(0) + xdata(:,:,2)*cos(0);
x2datarot(:,:,1)= xdata(:,:,1)*cos(0) - xdata(:,:,2)*sin(0);
x2datarot(:,:,2)= xdata(:,:,1)*sin(0) + xdata(:,:,2)*cos(0);

x10rot = x(2)*cos(0) - x(4)*sin(0);
y10rot = x(2)*sin(0) + x(4)*cos(0);

x20rot = x(8)*cos(0) - x(10)*sin(0);
y20rot = x(8)*sin(0) + x(10)*cos(0);

F = x(1)*exp(   -((x1datarot(:,:,1)-x10rot).^2/(2*x(3)^2) + (x1datarot(:,:,2)-y10rot).^2/(2*x(5)^2) )    )+ x(6)*exp(   -((x2datarot(:,:,1)-x20rot).^2/(2*x(8)^2) + (x2datarot(:,:,2)-y20rot).^2/(2*x(10)^2) )    );
% figure(3)
% alpha(0)
% imagesc(F)
% colormap('gray')
% figure(gcf)%bring current figure to front
% drawnow
% beep
% pause %Wait for keystroke