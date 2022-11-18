% plot scatter plot with mean +/- confidence interval
basal=[];
cLTP=[];

mean_basal=mean(basal);
mean_cLTP=mean(cLTP);
std_basal=std(basal);
std_cLTP=std(cLTP);
SEM_basal=std_basal/sqrt(length(basal));
SEM_cLTP=std_cLTP/sqrt(length(cLTP));
ci_basal=tinv(0.975,size(basal,1)-1);
yci_basal=bsxfun(@times,SEM_basal,ci_basal);
ci_cLTP=tinv(0.975,size(cLTP,1)-1);
yci_cLTP=bsxfun(@times,SEM_cLTP,ci_cLTP);

x_basal=1+0.1*randn(length(basal),1);
x_cLTP=2+0.1*randn(length(cLTP),1);
group_basal=ones(length(basal),1);
group_cLTP=2*ones(length(cLTP),1);
x=[x_basal;x_cLTP];
y=[basal;cLTP];
gs=[group_basal;group_cLTP];
%%
close all;
figure;
scat=gscatter(x,y,gs,"kr",'.');
scat.MarkerFace
hold on
plot([0.6; 1.4],mean_basal*ones(2,1),'k','LineWidth',2)
plot([0.7; 1.3],(yci_basal+mean_basal)*ones(2,1),'k','LineWidth',1)
plot([0.7; 1.3],(-yci_basal+mean_basal)*ones(2,1),'k','LineWidth',1)
plot([1,1],[yci_basal+mean_basal,-yci_basal+mean_basal],'k','LineWidth',1)

plot([1.6; 2.4],mean_cLTP*ones(2,1),'r','LineWidth',2)
plot([1.7; 2.3],(yci_cLTP+mean_cLTP)*ones(2,1),'r','LineWidth',1)
plot([1.7; 2.3],(-yci_cLTP+mean_cLTP)*ones(2,1),'r','LineWidth',1)
plot([2,2],[yci_cLTP+mean_cLTP,-yci_cLTP+mean_cLTP],'r','LineWidth',1)
ax=gca;

xlim([0.4 2.6]);xlabel('')
set(ax,'Xtick',[1,2],'XTickLabel',{'basal', 'cLTP'})
ylabel('# GluA1/synapse')
ylim([0 80])
legend('off')
