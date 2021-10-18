%Generates fig. 5 in the manuscript
%Requires fminsearchbnd.m
clear all
line='3013'
file=strcat('data/growthCurves_',line,'.mat')
load(file)
fun = @(x)distanceData(x(1),x(2),x(3),C);
x0 = [0.1,0.1,0.1];
fit = fminsearchbnd(fun,x0,[0 0 0])

hold on
plotFit(file,fit,1,'r','');
plotFit(file,fit,2,'g','');
plotFit(file,fit,3,'b','');
plotFit(file,fit,4,'k','');
plotFit(file,fit,5,'m','');
plotFit(file,fit,6,'y','');

%axis([0 50 0 0.03])
xlabel('time (hours)')
ylabel('normalised density')
set(gca,'FontSize',14)
legend('Allee model','data')

E=distanceData(fit(1),fit(2),fit(3),C)
n=200*6 %no of data points
AIC=2*3 + n*log(E)