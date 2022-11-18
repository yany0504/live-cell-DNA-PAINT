function F = Two2DGaussFunctionRot(x,xdata)
%% x = [Amp1, x01, wx1, y01, wy1, fi1, Amp2, x02, wx2, y02, wy2, fi2]
%[X,Y] = meshgrid(x,y) 
%  xdata(:,:,1) = X
%  xdata(:,:,2) = Y           
% Mrot = [cos(fi) -sin(fi); sin(fi) cos(fi)]
%%
x1datarot(:,:,1)= xdata(:,:,1)*cos(x(6)) - xdata(:,:,2)*sin(x(6));
x1datarot(:,:,2)= xdata(:,:,1)*sin(x(6)) + xdata(:,:,2)*cos(x(6));
x2datarot(:,:,1)= xdata(:,:,1)*cos(x(12)) - xdata(:,:,2)*sin(x(12));
x2datarot(:,:,2)= xdata(:,:,1)*sin(x(12)) + xdata(:,:,2)*cos(x(12));

x10rot = x(2)*cos(x(6)) - x(4)*sin(x(6));
y10rot = x(2)*sin(x(6)) + x(4)*cos(x(6));

x20rot = x(8)*cos(x(12)) - x(10)*sin(x(12));
y20rot = x(8)*sin(x(12)) + x(10)*cos(x(12));

F = x(1)*exp(   -((x1datarot(:,:,1)-x10rot).^2/(2*x(3)^2) + (x1datarot(:,:,2)-y10rot).^2/(2*x(5)^2) )    )+ x(7)*exp(   -((x2datarot(:,:,1)-x20rot).^2/(2*x(9)^2) + (x2datarot(:,:,2)-y20rot).^2/(2*x(11)^2) )    );
% figure(3)
% alpha(0)
% imagesc(F)
% colormap('gray')
% figure(gcf)%bring current figure to front
% drawnow
% beep
% pause %Wait for keystroke