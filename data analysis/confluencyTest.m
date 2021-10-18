%Calculates correlation between cell count and counfluency
clear all

load data/CO558.mat;
C3123 = CO(:,1:6);
cc1=mean(mean(C3123))
C3289 = CO(:,7:12);
cc2=mean(mean(C3289))

load data/CO560.mat;
C3013 = CO(:,1:6);
cc3=mean(mean(C3013))