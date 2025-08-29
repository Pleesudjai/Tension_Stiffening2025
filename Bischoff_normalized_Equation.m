clc; clear all; close all;

fc  = 48.8;
Eb  = 41000;
ft  = 2.6;
etu = 16000e-6;

% calculation
Em   = 4735*fc^0.5;
ecr  = ft/Em;
em   = linspace(ecr,etu,100);
beta = exp(-1100*(em-ecr)*Eb/200000);

strain = [0,em];
stress = [0,beta*ft];
beta   = [0,beta];


[strain;stress;beta]'

figure(1)
plot(strain,beta,'-r.');

figure(2)
plot(strain,stress,'-b.');