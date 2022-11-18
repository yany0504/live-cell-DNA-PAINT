%ANOVA YY
clear;close all; clc;
a = []; % manually put the data here
a(a == 0)= NaN;

[p,t,stats]=anova1(a); % anova
[c,m,h,cms]=multcompare(stats); % multiple comparison test
[c,m,h,cms]=multcompare(stats);
% number = size(a,2);
% for i = 1:number
%     for j=1:number
%         p(i,j)=anova1([a(:,i) a(:,j)],'');
%     end
% end
set(gca,'xtick',1:2, 'xticklabel',{'Synaptic','Juxtasynaptic'})
set(gca,'FontSize',20);
ylabel('Trajectory length (\mum, log)','FontSize',20);
title('Synaptic as 0~1 \mum')
hBox=figure(4);
set(hBox,'position',FigSize)