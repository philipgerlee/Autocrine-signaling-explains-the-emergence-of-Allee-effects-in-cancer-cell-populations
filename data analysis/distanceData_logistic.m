function f=distanceData_logistic(A,mu,C)
cut=200;
dt=0.25;
tmax=length(C(1,1:cut));

for i=1:6
    N(1)=C(i,1);
    for t=1:tmax-1
        N(t+1)=N(t)+dt*A*N(t)*(1-N(t))-dt*mu*N(t);
    end
    error(i)=norm(C(i,1:cut)-N);
end
f=mean(error);