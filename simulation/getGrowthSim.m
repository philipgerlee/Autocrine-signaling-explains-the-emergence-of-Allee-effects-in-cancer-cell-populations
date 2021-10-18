%Used for generating fig. 4
function f=getGrowthSim(mu,d)
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
d=factor*d; %diffusion coeff in cm^2/s
rho=factor*10*1e-5;
delta=factor*1*1e-5;
a=1*1e-5;
mu=mu*1e-5;
gamma=4*1.6e-10/h^2;

pc=0.01;
alpha=d*dt/h^2; %for the diff eq.


%IC
C=zeros(n);
g=zeros(n);

%Random plating
for i=1:n
    for j=1:n
        rr=rand;
        if rr<pc
            C(i,j)=1;
        else
            C(i,j)=0;
        end
    end
end


nc=zeros(1,tmax);
Nc=zeros(1,tmax);

Nc(1)=sum(sum(C))/n^2;

for t=1:tmax
    gh=(I-0.5*alpha*T)\(g*(I+0.5*alpha*T)+(0.5*dt*(rho*C-delta.*g)));
    g=((I+0.5*alpha*T)*gh+(0.5*dt*(rho*C-delta.*gh)))/(I-0.5*alpha*T);
        
    ind=find(C);
    for i=1:length(ind)
        ii=ind(i);
        wc=a*(1+g(ii));
        if rand<dt*wc %divide
            %long-range
            %xn=ceil(n*rand);
            %yn=ceil(n*rand);
            
            %short-range
            [x,y]=ind2sub([n n],ii);
            [xn,yn]=NN(x,y,n);
            if C(xn,yn)==0
                C(xn,yn)=1;
            end
        end
        if rand<dt*mu %die
            [x,y]=ind2sub([n n],ii);
            C(x,y)=0;
        end
        if rand<dt*gamma %move
            [x,y]=ind2sub([n n],ii);
            [xn,yn]=NN(x,y,n);
            if C(xn,yn)==0
                C(xn,yn)=1;
                C(x,y)=0;
            end
        end
    end
    nc(t)=sum(sum(C));
end


f=nc/n^2;






