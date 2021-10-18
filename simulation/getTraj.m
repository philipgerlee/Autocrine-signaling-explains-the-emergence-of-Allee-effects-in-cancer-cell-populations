%Required for generating fig. 3A
function [f,g]=getTraj(tmax,pc)
n=100; %number of cells
L=0.2; %size of the domain


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
d=factor*5e-11; %diffusion coeff in cm^2/s
rho=factor*10*1e-5; %production rate
delta=factor*1*1e-5; %decay rate
a=1*1e-5; %base proliferation rate
mu=3.75*1e-5; %death rate
gamma=0*4*2e-10/h^2; %jump rate

alpha=d*dt/h^2 %for the ADI-solver

%IC
C=zeros(n); %matrix for cells 
g=zeros(n); %matrix for GF

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
    %Solving the GF-equation
    gh=(I-0.5*alpha*T)\(g*(I+0.5*alpha*T)+(0.5*dt*(rho*C-delta.*g)));
    g=((I+0.5*alpha*T)*gh+(0.5*dt*(rho*C-delta.*gh)))/(I-0.5*alpha*T);
    
    ind=find(C);
    for i=1:length(ind)
        ii=ind(i);
        wc=a*(1+g(ii));
        if rand<dt*wc %divide
            %Long-range dispersal
            xn=ceil(n*rand);
            yn=ceil(n*rand);
            
            %Short-range dispersal
            %[x,y]=ind2sub([n n],ii);
            %[xn,yn]=NN(x,y,n);
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

S=1/(2*delta)+L/(4*sqrt(1*d*delta));
mus=a+a*rho*S/n

for t=1:tmax 
    Gc=a+a*rho*Nc(t)/delta + a*rho*(1-Nc(t))*S/n;
    Nc(t+1)=Nc(t)+dt*Gc*Nc(t)*(1-Nc(t))-dt*mu*Nc(t);
end

f=nc/n^2;
g=Nc;




