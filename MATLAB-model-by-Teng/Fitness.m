clear;
clc;

Col=linspecer(3);
p=2;%plasmid copy number;
K=[1 2 1+p];
c=0.1;
A=0:0.01:6;
n=3;


f1=(1-c)*K(1)^n./(K(1)^n+A.^n);
f2=(1-c)^2*K(2)^n./(K(2)^n+A.^n);
f3=(1-c)^(1+p)*K(3)^n./(K(3)^n+A.^n);

plot(A,f1,'color',Col(1,:),'linewidth',3);hold on;
plot(A,f2,'color',Col(2,:),'linewidth',3);hold on;
plot(A,f3,'color',Col(3,:),'linewidth',3);hold on;

legend({'f_1','f_2','f_3'},'FontAngle','italic');
legend boxoff;

set(gca,'fontsize',16);
xlabel('antibiotic concentration','fontsize',20);
ylabel('fitness','fontsize',20);
set(gcf,'position',[100 100 270 270]);

saveas(gcf,'Fitness.fig');
saveas(gcf,'Fitness.png');