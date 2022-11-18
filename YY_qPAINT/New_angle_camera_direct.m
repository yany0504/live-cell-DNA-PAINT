function New_angle_camera_direct(drift, pixel_size_ir,pixel_size_ccd,interp_factor)
%Corrects camera angle difference between IR and EMCCD cameras. 

%%
x1=drift(:,3);
y1=drift(:,4);
% z1=(drift(:,3)+0.497)/0.0023;
z1=(drift(:,5)-282.56633)/0.06768; % YY measured and 8/18/2020 
% z1=((drift(:,3)/4)+1.332)/0.000685; %matlab sigma=a*z+b => z=(sigma-b)/a
% xs = smooth(x1,0.1,'rloess');
% ys = smooth(y1,0.1,'rloess');
% zs = smooth(z1,0.1,'rloess');
xs = smooth(x1,0.1,'sgolay');
ys = smooth(y1,0.1,'sgolay');
zs = smooth(z1,0.1,'sgolay');
d =[xs,ys];



%from IR camera
% a1 = 681.8964*pixel_size_ir; %the first point
% b1 = 152.6751*pixel_size_ir;
% a2 = 638.2372*pixel_size_ir; %the last point
% b2 = 222.6751*pixel_size_ir;
% a1 = 897.000*pixel_size_ir; %YY 20200813
% b1 = 502.667*pixel_size_ir;
% a2 = 1016.333*pixel_size_ir; %YY 20200813
% b2 = 390.000*pixel_size_ir;
% a1 = 1018*pixel_size_ir; %YY 20200820
% b1 = 452*pixel_size_ir;
% a2 = 676*pixel_size_ir; %YY 20200820
% b2 = 469*pixel_size_ir;
a1 = 850*pixel_size_ir; %YY 20200820
b1 = 22*pixel_size_ir;
a2 = 763*pixel_size_ir; %YY 20200820
b2 = 26*pixel_size_ir;
L1 = sqrt((a2-a1)^2+(b2-b1)^2); %scale of vector

%from EMCCD camera
% c1 = 182.4220*pixel_size_ccd; %the first point
% d1 = 271.4277*pixel_size_ccd;
% c2 = 132.2694*pixel_size_ccd; %the last point
% d2 = 347.2374*pixel_size_ccd;
% c1 = 39.173*pixel_size_ccd; %YY 20200813
% d1 = 37.067*pixel_size_ccd;
% c2 = 44.533*pixel_size_ccd; %YY 20200813
% d2 = 32.240*pixel_size_ccd;
% c1 = 44.160*pixel_size_ccd; %YY 20200820
% d1 = 34.400*pixel_size_ccd;
% c2 = 29.493*pixel_size_ccd; %YY 20200820
% d2 = 35.200*pixel_size_ccd;
c1 = 364*pixel_size_ccd; %YY 20200820
d1 = 114*pixel_size_ccd;
c2 = 273*pixel_size_ccd; %YY 20200820
d2 = 120*pixel_size_ccd;
L2 = sqrt((c2-c1)^2+(d2-d1)^2); %scale of vector

%angle of two cameras
cos = ((a2-a1)*(c2-c1)+(b2-b1)*(d2-d1))/(L1*L2);
t = acosd(((a2-a1)*(c2-c1)+(b2-b1)*(d2-d1))/(L1*L2));

%% Plot box
figure
x = d(:,1);
y =d(:,2);
plot(x,y,'b*-');
shg


% create a 2D rotation matrix
rot = [cosd(t) sind(t);-sind(t) cosd(t)];

% Do Rotation
size=length(d);
B=zeros(size,2);
rot_pts=zeros(size,2);
zs_corr=zeros(size,1); 

for i=1:size;
  B = [d(i,1) d(i,2)];   
  Y = transpose(B);
  rot_pts(i,:) = rot*Y;
end

% from the original point
hold on
plot(rot_pts(:,1),rot_pts(:,2),'r*-');

for j=1:size;
    rot_pts_ori(j,1) = rot_pts(j,1)-rot_pts(1,1);
    rot_pts_ori(j,2) = rot_pts(j,2)-rot_pts(1,2);
    zs_corr(j) = zs(j) - zs(1);
end

table=[rot_pts_ori];

%interpolation the drift data
% interp_factor = 4 ; %interpolation factor
corr_x= interp(table(:,1),interp_factor);
corr_y= interp(table(:,2),interp_factor);
corr_z= interp(zs_corr,interp_factor);

drift_angle=[corr_x, corr_y, corr_z];  

assignin('caller','drift_angle',drift_angle);
assignin('base','drift_angle',drift_angle);


display('Correction completed.');
end



