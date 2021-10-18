%Based on data from single wells this script calculates growth curves
clear all
row={'A','B','C','D','E','F','G','H'};
col=1:6;
offset = 14
for i=1:length(col)
    D=zeros(8,370);
    for j=1:length(row)
        [i j]
        load(strcat('data/560_',row{j},num2str(col(i)),'_confluency.mat'),'confluency');
        load(strcat('data/560_',row{j},num2str(col(i)),'_times.mat'),'times');
        
        up=times(end)
        nop=round(times(end)/15) 
        Tinp=linspace(0,up,nop)
        
        conf=interp1(times,confluency,Tinp)
        conf=conf(offset:end)
        
        D(j,1:length(conf))=conf;
    end
    C(i,:)=mean(D);
end
plot(Tinp(offset:end),C','LineWidth',2)
xlabel('time')
ylabel('normalised density')
