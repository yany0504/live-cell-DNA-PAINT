close all;
%figure;
output=scatplot(Diff_tot,Traj_tot,'N',12);
xx=output.xi;
yy=output.yi;
zz=output.zi/max(max(output.zi));
gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(GMModel_tot,[x0 y0]),x,y);
z_value=gmPDF(xx,yy);

figure; 
ax1=subplot(1,3,1);
surf(xx,yy,zz);
xlabel('log diff. coeff, (log(\mum^2/s))');
ylabel('log traj. range, (log(\mum))');
zlabel('PDF');
title('Raw data');
%figure;
% S=fsurf(gmPDF,[min(min(xx)), max(max(xx)), min(min(yy)),max(max(yy))]);
residual=zz-z_value;
ax2=subplot(1,3,2);surf(xx,yy,z_value);
xlabel('log diff. coeff, (log(\mum^2/s))');
ylabel('log traj. range, (log(\mum))');
zlabel('PDF');
title('GM model');
ax3=subplot(1,3,3);surf(xx,yy,residual);
xlabel('log diff. coeff, (log(\mum^2/s))');
ylabel('log traj. range, (log(\mum))');
zlabel('PDF');
title('Residual');
Link = linkprop([ax1, ax2, ax3],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
set(gca,'Zlim',[-0.5,1.2]);
setappdata(gcf, 'StoreTheLink', Link);
set(gcf,'Position',[100 100 1500 500]);