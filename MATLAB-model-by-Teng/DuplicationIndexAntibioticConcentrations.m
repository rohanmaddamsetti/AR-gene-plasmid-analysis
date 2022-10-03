clear;
clc;

global c A K n p D eta;
p=2;%plasmid copy number;
K=[1 2 1+p];
cs=0.05:0.05:0.25;
As=[0.25:0.01:1.2];
n=3;
D=0.1;
eta=2*10^(-4);
Col=linspecer(length(cs));

initial=[1 0 0];
timespan=0:200;
for i=1:length(cs)
    
    c=cs(i);
    DuplicationIndex=0*ones(length(As),0);
    for j=1:length(As)
        j
        A=As(j);
        [t,y]=ode45(@GeneDuplication,timespan,initial);
        Duplication(j)=(y(end,2)+y(end,3))./sum(y(end,:),2);
    end
    plot(As,Duplication,'color',Col(i,:),'linewidth',3);hold on;
end
axis([min(As) max(As) 0 1]);
set(gca,'fontsize',16);
xlabel('antibiotic concentration','fontsize',20);
ylabel('duplication index','fontsize',20);
set(gcf,'position',[100 100 270 270]);

saveas(gcf,'DuplicationIndexAntibioticConcentration.fig');
saveas(gcf,'DuplicationIndexAntibioticConcentration.png');

function dydt=GeneDuplication(t,y)
global c A K n p D eta;
f1=(1-c)*K(1)^n/(K(1)^n+A^n);
f2=(1-c)^2*K(2)^n/(K(2)^n+A^n);
f3=(1-c)^(1+p)*K(3)^n/(K(3)^n+A^n);
yt=y(1)+y(2)+y(3);
dydt=[f1*y(1)*(1-yt)-D*y(1)-2*eta*y(1);
    f2*y(2)*(1-yt)-D*y(2)+eta*y(1);
    f3*y(3)*(1-yt)-D*y(3)+eta*y(1)];
end
