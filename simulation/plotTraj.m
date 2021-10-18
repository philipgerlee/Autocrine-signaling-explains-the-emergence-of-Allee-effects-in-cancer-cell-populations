%Generates fig. 3A in the manuscript
clear all
tmax=1e4;
dt=1e2;
pc=[0.05 0.3 0.8];

for i=1:length(pc)
    i
    [ns na]=getTraj(tmax,pc(i));
    S(:,i)=ns;
    A(:,i)=na;
end

T=dt*(1:tmax);
hold on
for i=1:length(pc)
    plot(T,S(:,i),'k','LineWidth',2)
    plot(T,A(1:end-1,i),'k:','LineWidth',2)
end
ylabel('normalised density')
xlabel('time (seconds)')
ylim([0 1])
set(gca,'FontSize',14)
legend('IB-model','ODE-model')
