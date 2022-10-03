clear;
clc;

global c A K n p D eta;
p=2;%plasmid copy number;
K=[1 2 1+p];
c=0.1;
A=2;
n=3;
D=0.1;
etas=[0 2*10^(-7) 2*10^(-6) 2*10^(-5) 2*10^(-4)];
Col=linspecer(length(etas));

initial=[1 0 0];
timespan=0:200;
for i=1:length(etas)
    eta=etas(i);
[t,y]=ode45(@GeneDuplication,timespan,initial);

plot(t,(y(:,2)+y(:,3))./sum(y,2),'color',Col(i,:),'linewidth',3);hold on;
end
set(gca,'fontsize',16);
xlabel('time','fontsize',20);
ylabel('duplication index','fontsize',20);
set(gcf,'position',[100 100 270 270]);

saveas(gcf,'DuplicationIndexTimeCourse.fig');
saveas(gcf,'DuplicationIndexTimeCourse.png');

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
