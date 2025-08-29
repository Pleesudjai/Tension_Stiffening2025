function mtxST = Normal_Dist_Strength(ckPat,avg,std,minST,maxST,stateNo,nNode,endL,clearL,minSpac,h)

% Crack patterns
% 1) binomial: such that crack start at 50, [25,75], [12.5,37.5,62.5,87.5] ....
% 2) weak at center, such that crack start at 50, [50-minSpac,50+minSpac], [50-2minSpac,50+2minSpac],...
% 3) from left to right, such that crack start from 0,minSpac,2minSpac,3minSpac...
% 4) random: crack start at any random locations
mtxST = linspace(maxST,maxST,nNode);
if ckPat==1 || ckPat==2
    i      = 1;
    sep(1) = 1;
    if ckPat==1
        % 1 deterministic crack, symmetry strength zig zac
        x(1)    = endL+clearL/2;
        ndiv(1) = 2;
        subL    = clearL;
        while subL > minSpac
            i    = i+1;
            ndiv(i) = 2^i;
            subL = clearL/ndiv(i);
            xi   = endL + linspace(subL,clearL-subL,ndiv(i-1));
            x    = [x,xi];
            sep(i) = length(x);
        end
    elseif ckPat==2
        % 2 deterministic crack, symmetrical strength distribution with the weakest at the middle and
        % gradually increase toward the end grip
        xm   = endL + clearL/2;
        x(1) = xm;
        di   = 0;
        dMax = clearL/2
        while di < dMax
            i     = i+1;
            di    = (i-1)*minSpac;
            x     = [x,xm-di,xm+di];
            sep(i)= length(x);
        end
    end
    nGroup = i;
    node   = round(x/h)+1;
    beg    = 1;
    prob   = sort(normrnd(avg,std,length(x),1)');
    prob(1) = prob(1)*0.9999;    % lower the strength at the center such that the first crack appear at the center
    for i=1:nGroup;
        mtxST(node(beg:sep(i))) = mean(prob(beg:sep(i)));
        %mtxST(node(beg:sep(i))) = prob(beg:sep(i));
        beg  = sep(i)+1;
    end
else
    nPoint = round(clearL/minSpac)+1;
    randn('state', stateNo);
    if ckPat==3 % 3 deterministic crack, min strength from left to right
        prob   = sort(normrnd(avg,std,nPoint,1)');
    elseif ckPat==4 % 4 random crack
        prob   = normrnd(avg,std,nPoint,1)';
    elseif ckPat==5  % special purpose
        prob = normrnd(avg,std,nPoint,1)';
        minProb = min(prob);
        mid  = round(length(prob)/2);
        prob(mid)     = minProb*0.9995;  % lower the strength at center such that the first crack appear at center
        prob([1,end]) = minProb*0.9999; % lower the strength at each end of the clear span such that the second crack appear at each end clear span
    end
    mtxST = linspace(maxST,maxST,nNode);
    x     = linspace(endL,endL+clearL,nPoint);
    node  = round(x/h)+1;
    nbeg  = round(endL/minSpac)+1;
    for i=1:nPoint
        if prob(i) < minST
            prob(i) = minST;
        elseif prob(i) > maxST
            prob(i) = maxST;
        end
        mtxST(node(i)) = prob(i);
    end  
end