%Required for fig. 5
%Plots data and model fit
function f=plotFit(file,P,column,color,dash)

load(file)
dt=0.25; %in hours
A=P(1);
B=P(2);    
mu=P(3);
cut=200;

tmax=length(C(1,1:cut));
T=dt*(1:length(C(1,1:cut)));
time=dt*(1:tmax);

N(1)=C(column,1);
for t=1:tmax-1 
    N(t+1)=N(t)+dt*(A+B*N(t))*N(t)*(1-N(t))-dt*mu*N(t);
end

hold on
plot(time,N,strcat(color,dash,'-'),'LineWidth',2)
plot(T,C(column,1:cut),strcat(color,'o'),'LineWidth',2)
f=[];