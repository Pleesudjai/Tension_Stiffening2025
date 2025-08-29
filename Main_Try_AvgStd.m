clear all;
clc;

% input
nSample = 5;
avg     = 3;
std     = 0.4;
clearL  = 200;
minSpac = 2;

% calculation
nPoint = round(clearL/minSpac)+1;
for i=1:nSample
    randn('state', i);
    prob   = normrnd(avg,std,nPoint,1);
    minProb(i) = min(prob);
    maxProb(i) = max(prob);
end
minProb
minfcr = mean(minProb)
maxProb
maxfcr = mean(maxProb)