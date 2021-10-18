%This script generates fig. 2 in the manuscript
clear all

%parameters
L=.2;
factor=100;
d=factor*5e-11; 
rho=factor*10*1e-5;
delta=factor*1*1e-5;
a=1*1e-5;
n=100;

S=1/(2*delta)+L/(4*sqrt(1*d*delta));

N=linspace(0,1,1000);

%%%%%
mu=3.286*1e-5;
Gamma=a+a*rho*N/delta + a*rho*(1-N)*S/n;
G=Gamma.*N.*(1-N)-mu*N;
hold on
plot(N,G./N,'LineWidth',2) %growth rate
%plot(N,G,'b','LineWidth',2) %per-capita growth rate

%%%%%%
mu=3.85*1e-5;
S=1/(2*delta)+L/(4*sqrt(1*d*delta));
D=(4*mu*delta*rho*n*(delta*S-n)+a*(rho*n-delta*(n+2*rho*S))^2)/(4*a*(n-delta*S)^2*rho^2);
E=a+a*rho*S/n-mu;

Gamma=a+a*rho*N/delta + a*rho*(1-N)*S/n;
G=Gamma.*N.*(1-N)-mu*N;
hold on
plot(N,G./N,'r','LineWidth',2)
%plot(N,G,'r','LineWidth',2)

%%%%%%
mu=0*1e-5;
S=1/(2*delta)+L/(4*sqrt(1*d*delta));


Gamma=a+a*rho*N/delta + a*rho*(1-N)*S/n;
G=Gamma.*N.*(1-N)-mu*N;
plot(N,G./N,'g','LineWidth',2)
%plot(N,G,'g','LineWidth',2)

plot(N,zeros(1,length(N)),'k--')

xlabel('population density')
ylabel('per capita growth rate f(n)/n')
%ylabel('growth rate f(n)')
legend('\mu = \mu_c','\mu = 1.15 \mu_c','\mu=0')
xlim([0 1])
%ylim([-0.1 0.3]*1e-5)
set(gca,'FontSize',28)
set(gca,'FontSize',14)

D=(4*mu*delta*rho*n*(delta*S-n)+a*(rho*n-delta*(n+2*rho*S))^2)/(4*a*(n-delta*S)^2*rho^2)
