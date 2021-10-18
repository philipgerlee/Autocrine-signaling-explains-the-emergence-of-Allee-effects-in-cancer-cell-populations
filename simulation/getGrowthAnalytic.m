%Used for generating fig. 4
function g=getGrowthAnalytic(mu,d)
L=0.2; %size of domain
n=100; %number of cells
tmax=1e4; %duration of simulation


%Initialising the Laplacian matrix
T=diag(-2*ones(1,n)) + diag(ones(n-1,1),1) + diag(ones(n-1,1),-1);
T(1,1)=-1;
T(n,n)=-1;
I=eye(n);

%Numerics
h=L/n; %cell size in cm
dt=1e2; %time step in seconds

factor=100;
%Parameters
d=factor*d;
rho=factor*10*1e-5;
delta=factor*1*1e-5;
a=1*1e-5;
mu=mu*1e-5;
gamma=4*1e-8/h^2;

pc=0.01;
alpha=d*dt/h^2; %for the diff eq.

Nc=zeros(1,tmax);

Nc(1)=pc;

K=1/(2*delta)+L/(4*sqrt(1*d*delta));

for t=1:tmax 
    Gc=a+a*rho*Nc(t)/delta + a*rho*(1-Nc(t))*K/n;
    Nc(t+1)=Nc(t)+dt*Gc*Nc(t)*(1-Nc(t))-dt*mu*Nc(t);
end
g=Nc;







