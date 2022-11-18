% making plot of diffusion coefficient as a function of distance from
% Homer1
% figure;
% Output=scatplot(Dist_tot,Diff_tot);
% xlabel('Distance from NN-Homer1 (\mum)');
% ylabel('Diffusion coefficient (log, \mum^2/s)');
% title('Total')
FigSize=[100, 100, 800, 600];
figure;
Output=scatplot(Diff_tot,Dist_tot);
ylabel('Distance from NN-Homer1 (\mum)','FontSize',25);
xlabel('Diffusion coefficient (log, \mum^2/s)','FontSize',25);
% title('Total')
set(gca,'Color',[0.24 0.15 0.6],'FontSize',20);
xlim([-5 -0.5]);
ylim([-0.05 2.05]);
syn_immo=Diff_tot(Diff_tot(:,1)<-2.5&Dist_tot<0.4);
set(gcf, 'Position', FigSize);
% figure;
% Output=scatplot(Dist_syn,Diff_syn);
% xlabel('Distance from NN-Homer1 (\mum)');
% ylabel('Diffusion coefficient (log, \mum^2/s)');
% title('Synaptic')
% figure;
% Output=scatplot(Diff_syn,Dist_syn);
% ylabel('Distance from NN-Homer1 (\mum)');
% xlabel('Diffusion coefficient (log, \mum^2/s)');
% title('Synaptic')
% set(gca,'Color',[0.24 0.15 0.6]);
% xlim([-5.5 -0.5]);
% figure;
% Output=scatplot(Dist_juxta,Diff_juxta);
% xlabel('Distance from NN-Homer1 (\mum)');
% ylabel('Diffusion coefficient (log, \mum^2/s)');
% title('Juxta-synaptic')
% figure;
% Output=scatplot(Diff_juxta,Dist_juxta);
% ylabel('Distance from NN-Homer1 (\mum)');
% xlabel('Diffusion coefficient (log, \mum^2/s)');
% title('Juxta-synaptic')
% set(gca,'Color',[0.24 0.15 0.6]);
% xlim([-5.5 -0.5]);


dist_diff=[Dist_tot, Diff_tot];



diff_100=dist_diff(dist_diff(:,1)<=0.1,2);
diff_200=dist_diff(dist_diff(:,1)<=0.2&dist_diff(:,1)>0.1,2);
diff_300=dist_diff(dist_diff(:,1)<=0.3&dist_diff(:,1)>0.2,2);
diff_400=dist_diff(dist_diff(:,1)<=0.4&dist_diff(:,1)>0.3,2);
diff_500=dist_diff(dist_diff(:,1)<=0.5&dist_diff(:,1)>0.4,2);
diff_600=dist_diff(dist_diff(:,1)<=0.6&dist_diff(:,1)>0.5,2);
diff_700=dist_diff(dist_diff(:,1)<=0.7&dist_diff(:,1)>0.6,2);
diff_800=dist_diff(dist_diff(:,1)<=0.8&dist_diff(:,1)>0.7,2);
diff_900=dist_diff(dist_diff(:,1)<=0.9&dist_diff(:,1)>0.8,2);
diff_1000=dist_diff(dist_diff(:,1)<=1.0&dist_diff(:,1)>0.9,2);
diff_1100=dist_diff(dist_diff(:,1)<=1.1&dist_diff(:,1)>1.0,2);
diff_1200=dist_diff(dist_diff(:,1)<=1.2&dist_diff(:,1)>1.1,2);
diff_1300=dist_diff(dist_diff(:,1)<=1.3&dist_diff(:,1)>1.2,2);
diff_1400=dist_diff(dist_diff(:,1)<=1.4&dist_diff(:,1)>1.3,2);
diff_1500=dist_diff(dist_diff(:,1)<=1.5&dist_diff(:,1)>1.4,2);
diff_1600=dist_diff(dist_diff(:,1)<=1.6&dist_diff(:,1)>1.5,2);
diff_1700=dist_diff(dist_diff(:,1)<=1.7&dist_diff(:,1)>1.6,2);
diff_1800=dist_diff(dist_diff(:,1)<=1.8&dist_diff(:,1)>1.7,2);
diff_1900=dist_diff(dist_diff(:,1)<=1.9&dist_diff(:,1)>1.8,2);
diff_2000=dist_diff(dist_diff(:,1)<=2.0&dist_diff(:,1)>1.9,2);


%% calculate immobile/total
criteria=-2.5;
x_dist=[0.1:0.1:2];
y_immo_frac=[length(diff_100(diff_100<criteria))/length(diff_100),length(diff_200(diff_200<criteria))/length(diff_200),length(diff_300(diff_300<criteria))/length(diff_300),length(diff_400(diff_400<criteria))/length(diff_400),...
    length(diff_500(diff_500<criteria))/length(diff_500),length(diff_600(diff_600<criteria))/length(diff_600),length(diff_700(diff_700<criteria))/length(diff_700),length(diff_800(diff_800<criteria))/length(diff_800),...
    length(diff_900(diff_900<criteria))/length(diff_900),length(diff_1000(diff_1000<criteria))/length(diff_1000),length(diff_1100(diff_1100<criteria))/length(diff_1100),length(diff_1200(diff_1200<criteria))/length(diff_1200),...
    length(diff_1300(diff_1300<criteria))/length(diff_1300),length(diff_1400(diff_1400<criteria))/length(diff_1400),length(diff_1500(diff_1500<criteria))/length(diff_1500),length(diff_1600(diff_1600<criteria))/length(diff_1600),...
    length(diff_1700(diff_1700<criteria))/length(diff_1700),length(diff_1800(diff_1800<criteria))/length(diff_1800),length(diff_1900(diff_1900<criteria))/length(diff_1900),length(diff_2000(diff_2000<criteria))/length(diff_2000)];
figure;

plot(x_dist,y_immo_frac);
set(gca,'FontSize',20);
xlabel('Diffusion coefficient (log, \mum^2/s)','FontSize',25);
ylabel('Fraction immobile GluA1-AMPAR','FontSize',25);


%% plot every window

figure;
h100=histogram(diff_100,'BinWidth',0.2,'Normalization','Probability');
xlim([-6 -0.5]);
xlabel('log(diffusion coefficient) (\mum^2/s)');
ylabel('Frequency');
title('Distance 0~100 nm')

figure;
h200=histogram(diff_200,'BinWidth',0.2,'Normalization','Probability');
xlim([-6 -0.5]);
xlabel('log(diffusion coefficient) (\mum^2/s)');
ylabel('Frequency');
title('Distance 100~200 nm')

figure;
h300=histogram(diff_300,'BinWidth',0.2,'Normalization','Probability');
xlim([-6 -0.5]);
xlabel('log(diffusion coefficient) (\mum^2/s)');
ylabel('Frequency');
title('Distance 200~300 nm')

figure;
h400=histogram(diff_400,'BinWidth',0.2,'Normalization','Probability');
xlim([-6 -0.5]);
xlabel('log(diffusion coefficient) (\mum^2/s)');
ylabel('Frequency');
title('Distance 300~400 nm')

figure;
h500=histogram(diff_500,'BinWidth',0.2,'Normalization','Probability');
xlim([-6 -0.5]);
xlabel('log(diffusion coefficient) (\mum^2/s)');
ylabel('Frequency');
title('Distance 400~500 nm')

figure;
h600=histogram(diff_600,'BinWidth',0.2,'Normalization','Probability');
xlim([-6 -0.5]);
xlabel('log(diffusion coefficient) (\mum^2/s)');
ylabel('Frequency');
title('Distance 500~600 nm')

figure;
h700=histogram(diff_700,'BinWidth',0.2,'Normalization','Probability');
xlim([-6 -0.5]);
xlabel('log(diffusion coefficient) (\mum^2/s)');
ylabel('Frequency');
title('Distance 600~700 nm')

figure;
h800=histogram(diff_800,'BinWidth',0.2,'Normalization','Probability');
xlim([-6 -0.5]);
xlabel('log(diffusion coefficient) (\mum^2/s)');
ylabel('Frequency');
title('Distance 700~800 nm')

figure;
h900=histogram(diff_900,'BinWidth',0.2,'Normalization','Probability');
xlim([-6 -0.5]);
xlabel('log(diffusion coefficient) (\mum^2/s)');
ylabel('Frequency');
title('Distance 800~900 nm')

figure;
h1000=histogram(diff_1000,'BinWidth',0.2,'Normalization','Probability');
xlim([-6 -0.5]);
xlabel('log(diffusion coefficient) (\mum^2/s)');
ylabel('Frequency');
title('Distance 900~1000 nm')

figure;
h1100=histogram(diff_1100,'BinWidth',0.2,'Normalization','Probability');
xlim([-6 -0.5]);
xlabel('log(diffusion coefficient) (\mum^2/s)');
ylabel('Frequency');
title('Distance 1000~1100 nm')

figure;
h1200=histogram(diff_1200,'BinWidth',0.2,'Normalization','Probability');
xlim([-6 -0.5]);
xlabel('log(diffusion coefficient) (\mum^2/s)');
ylabel('Frequency');
title('Distance 1100~1200 nm')

figure;
h1300=histogram(diff_1300,'BinWidth',0.2,'Normalization','Probability');
xlim([-6 -0.5]);
xlabel('log(diffusion coefficient) (\mum^2/s)');
ylabel('Frequency');
title('Distance 1200~1300 nm')

figure;
h1400=histogram(diff_1400,'BinWidth',0.2,'Normalization','Probability');
xlim([-6 -0.5]);
xlabel('log(diffusion coefficient) (\mum^2/s)');
ylabel('Frequency');
title('Distance 1300~1400 nm')

figure;
h1500=histogram(diff_1500,'BinWidth',0.2,'Normalization','Probability');
xlim([-6 -0.5]);
xlabel('log(diffusion coefficient) (\mum^2/s)');
ylabel('Frequency');
title('Distance 1400~1500 nm')

figure;
h1600=histogram(diff_1600,'BinWidth',0.2,'Normalization','Probability');
xlim([-6 -0.5]);
xlabel('log(diffusion coefficient) (\mum^2/s)');
ylabel('Frequency');
title('Distance 1500~1600 nm')

figure;
h1700=histogram(diff_1700,'BinWidth',0.2,'Normalization','Probability');
xlim([-6 -0.5]);
xlabel('log(diffusion coefficient) (\mum^2/s)');
ylabel('Frequency');
title('Distance 1600~1700 nm')

figure;
h1800=histogram(diff_1800,'BinWidth',0.2,'Normalization','Probability');
xlim([-6 -0.5]);
xlabel('log(diffusion coefficient) (\mum^2/s)');
ylabel('Frequency');
title('Distance 1700~1800 nm')

figure;
h1900=histogram(diff_1900,'BinWidth',0.2,'Normalization','Probability');
xlim([-6 -0.5]);
xlabel('log(diffusion coefficient) (\mum^2/s)');
ylabel('Frequency');
title('Distance 1800~1900 nm')

figure;
h2000=histogram(diff_2000,'BinWidth',0.2,'Normalization','Probability');
xlim([-6 -0.5]);
xlabel('log(diffusion coefficient) (\mum^2/s)');
ylabel('Frequency');
title('Distance 1900~2000 nm')

