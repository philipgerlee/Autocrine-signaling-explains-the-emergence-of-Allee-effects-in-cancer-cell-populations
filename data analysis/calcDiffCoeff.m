%Calculates the diffusion coefficient from data
clear all
load data/558_MSD_slope.mat

D3123 = mean(MSD_slope(1:6));
D3289 = mean(MSD_slope(7:12));

load data/560_MSD_slope.mat

D3013 = mean(MSD_slope(1:6));

(D3123/4)*2.8e-12
(D3289/4)*2.8e-12
(D3013/4)*2.8e-12
