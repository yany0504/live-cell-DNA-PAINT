Names={'5xR1', 'R1'};
kon=[3.99*10^6 3*10^5];
stdkon=[2.59*10^6 2*10^5];
koff=[1.82033 10.05];
stdkoff=[0.570059 1.12];


figure;
hold on
yyaxis left
ylabel('k_{on} (M^{-1}s^{-1})','FontSize',25,'fontweight','bold')
bar([1 3],kon,0.3,'k');
er1=errorbar([1 3],kon,stdkon,'Marker','o','Color','k','LineStyle','none');
set(gca,'ycolor','k','FontSize',20);
yyaxis right
ylabel('k_{off} (s^{-1})','FontSize',25,'fontweight','bold');
bar([2 4],koff,0.3,'r');
er2=errorbar([2 4],koff,stdkoff,'Color','r','Marker','o','LineStyle','none');
set(gca,'ycolor','r','FontSize',20);
set(gca,'xtick',[1.5 3.5],'xticklabel',Names,'xcolor','k','FontSize',20,'fontweight','bold');

xlim([0.5 4.5]);
