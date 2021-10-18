%Generates fig. 5 in the manuscript
%Requires fminsearchbnd.m
clear all
line='3289'
file=strcat('data/growthCurves_',line,'.mat')
load(file)
fun = @(x)distanceData_logistic(x(1),x(2),C);
x0 = [0.2,0.2];
fit = fminsearchbnd(fun,x0,[0 0])


hold on
plotFit(file,[fit(1) 0 fit(2)],1,'r','-');
plotFit(file,[fit(1) 0 fit(2)],2,'g','-');
plotFit(file,[fit(1) 0 fit(2)],3,'b','-');
plotFit(file,[fit(1) 0 fit(2)],4,'k','-');
plotFit(file,[fit(1) 0 fit(2)],5,'m','-');
plotFit(file,[fit(1) 0 fit(2)],6,'y','-');

%axis([0 50 0 0.6])
xlabel('time (hours)')
ylabel('normalised density')
set(gca,'FontSize',14)
legend('logistic model','data')

E=distanceData_logistic(fit(1),fit(2),C)
n=200*6 %no of data points
AIC=2*2 + n*log(E)
