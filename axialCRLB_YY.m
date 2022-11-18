%% Calculating axial CRLB based on Rieger, B., & Stallinga, S. (2014). The Lateral and Axial Localization Uncertainty in Super?Resolution Light Microscopy. ChemPhysChem, 15(4), 664-670.
% You need N: signal photons, B: BG photons, sigma_0: spot size, which can
% be got from ThunderSTORM (qPALM maybe also possible)
% a: pixel size, l: 1/2 of the distance between vertical and horizontal sigma, 
% d: the focal depth (T-LSFM setup with 60x NA1.30 SiI objective, 0.58 um)
% output: deltaz=(l^2+d^2)/2lsqrt(N)*sqrt(1+8*tau+sqrt((9*tau)/(1+4*tau)))
% where tau = 2*pi*B*(sigma_0^2*(1+l^2/d^2)+a^2/12)/(N*(a^2))
%% Loading file
path = 'E:\Data\DNA-PAINT\200821_DNA-PAINT_repeat\2\';
if ~exist('fileName','var')|| isempty(fileName)
    [userfilein, userdirin]=uigetfile({
         '*.csv','Data file (*.csv)';...
        '*.*','All Files (*.*)'},'Select the ThunderSTORM result file to process',...
        path);
    fileName=fullfile(userdirin,userfilein);
else
    if ~exist(fileName,'file')
        fprintf('File not found: %s\n',fileName);
        return;
    end
end
data=xlsread(fileName);
%% analysis
N = 2000;%data(:,8);
b = 2;%data(:,10);
z=-300:300;
% sigma_x = 1.901933e-06*(z+472.006).^2 +0.000000e+00*(z+472.006).^3 +1.39438;
% sigma_y = 7.990730e-07*(z-670.626).^2 +0.000000e+00*(z-670.626).^3 +1.38260;
% sigma_x = 5.160056e-06*(z-334.832).^2 +1.16699;
% sigma_y = 3.463602e-06*(z+237.747).^2 +1.56630;
sigma_x = 6.690634e-06*(z-139.354).^2 +0.000000e+00*(z-139.354).^3 +1.23400;
sigma_y = 6.399794e-06*(z+170.845).^2 +0.000000e+00*(z+170.845).^3 +1.16169;

sigma_x=sigma_x*122;
sigma_y=sigma_y*122;
sigma_0=(sigma_x+sigma_y)./2;
a=107; %pixel size in nm
%l=sqrt(294.8066723164101*200.38705676842548); % 1/2 of the distance
l=150;
% between vertical and horizontal sigma in calibration curve. from the
%calibration function l = sqrt(c1*c2)
d=600*1.33/(1.5)^2;

tau = 2*pi.*(b.^2).*(sigma_x.*sigma_y.*(1+(l^2/d^2))+a^2/12)./(N.*(a^2));
F2=(4.*l^2.*z.^2)./((l^2+d^2+z.^2).^2);
deltasigma_x=(2.*sigma_x.^2+a^2/12)./N.*(1+8.*tau+sqrt((9.*tau)./(1+4.*tau)));
deltasigma_y=(2.*sigma_y.^2+a^2/12)./N.*(1+8.*tau+sqrt((9.*tau)./(1+4.*tau)));

delta_F2=(1-F2).*((deltasigma_x./sigma_x).^2+(deltasigma_y./sigma_y).^2);

deltaz=delta_F2.*((l^2+d^2+z.^2).^4./(4*l^2*(l^2+d^2-z.^2).^2));
deltaz=sqrt(deltaz);

deltax=sqrt(((sigma_0.^2+(a^2/12))./N).*(1+4.*tau+sqrt((2.*tau)./(1+tau))));

deltasigma_x=sqrt(deltasigma_x);
deltasigma_y=sqrt(deltasigma_y);





figure;
hxy=histogram(deltax);
xlabel('Lateral CRLB (nm)');
ylabel('Count');
figure;
hz=histogram(deltaz);
xlabel('Axial CRLB (nm)');
ylabel('Count');
median(deltax)
mean(deltax)
median(deltaz)
mean(deltaz)
