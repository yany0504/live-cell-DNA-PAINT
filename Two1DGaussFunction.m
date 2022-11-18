function F = Two1DGaussFunction(x,xdata)
%% x = [Amp1, x01, wx1,Amp2, x02, wx2]
%  xdata = X

%%


F = x(1)*exp(-((xdata-x(2)).^2/(2*x(3)^2)))+ x(4)*exp(-((xdata-x(5)).^2/(2*x(6)^2)));
% figure(3)
% alpha(0)
% imagesc(F)
% colormap('gray')
% figure(gcf)%bring current figure to front
% drawnow
% beep
% pause %Wait for keystroke