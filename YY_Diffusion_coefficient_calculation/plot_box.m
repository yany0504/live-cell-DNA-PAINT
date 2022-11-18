d1=Dif_Inv605_CTX_2345_Inv605; 
d2=Dif_NN620_0304;
dd=[d1;d2];
group=[ones(1,length(d1)), 2*ones(1,length(d2))]
boxplot(dd, group)