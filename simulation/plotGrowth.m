%For generating fig. 3B and 4
clear all
mu_s=linspace(2,4,15); %range of death rates
mu_a=linspace(2,4,200);

for i=1:length(mu_s)
    i
    ns=getGrowthSim(mu_s(i),5e-11);
    S(i)=ns(end);
end

for i=1:length(mu_a)
    i
    na=getGrowthAnalytic(mu_a(i),5e-11);
    A(i)=na(end);
end

factor=100;
L=.2;
d=factor*5e-11; %diffusion coeff in cm^2/s
rho=factor*10*1e-5; %production rate
delta=factor*1*1e-5; %decay rate
a=1*1e-5; %baseline proliferation rate
K=1/(2*delta)+L/(4*sqrt(d*delta));
n=100; %number of cells

mus=a+a*rho*K/n;
hold on
plot(mu_s,S,'ko','LineWidth',2)
plot(mu_a,A,'k','LineWidth',2)
plot(mus*ones(1,100)/1e-5,linspace(0,0.7,100),'k--','LineWidth',2)
xlabel('\mu (\times 10^5)')
ylabel('density after 11 days')
legend('Simulation','Analytical')
set(gca,'FontSize',14)