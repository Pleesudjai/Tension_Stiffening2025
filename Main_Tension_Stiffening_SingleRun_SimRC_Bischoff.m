close all;  clear all;  clc;    

% ------------------------------------
% Part A: Input data
% ------------------------------------
progVersion  = 'Main_Tension_Stiffening_SingleRun_SimRC_Bischoff 08/12/2007';   % program version
fname        = 'RND_STL_bond13fcr26rho2_Gen07.dat';                             % file name to report the output
fname_ckResp = 'RND_STL_bond13fcr26rho2_Gen07_ckResp.dat'                       % file name to report the responses at each crack formation
 
scaling   = 1.0;                                % scaling for yITF
avgST     = 2.6;                                % average matrix strength
totLoad   = 2*5e4;                              % total applied tensile load
rho       = 0.02;                               % reinforcement ratio
stateNo   = 7;                                  % state number to generate random number; need for ckPat=3,4,5
stdST     = 0.02;                               % standard deviation of matrix strength
 
% specimen geometry
endL    = 100;                                  % embedment length at the end grips
measL   = 900;                                  % embedment length for measurement in the middle
L       = endL + measL + endL;                  % total embedment length
 
% load and discretization
%totLoad = 3*5e4;                               % total applied tensile load
nNode   = 1101;                                 % number of nodes to approximate total embedment length, L
nInc    = 301;                                  % number of load increments
nRecPt  = 10;                                   % number of nodes used in recording distribution value at each increment, nRecPt <= nNode
dashType  = strvcat('--r','--m','--y','--k','--g','--c','--b');     % different colors used in ploting for selected increment, make sure that number of color > number of crack formation
solidType = strvcat('-r','-m','-y','-k','-g','-c','-b');            % different colors used in ploting for selected increment, make sure that number of color > number of crack formation
 
% interface shear strength model
xITF   = [0,  0.01, 0.11,  0.20,  0.42,  0.75,  1.20,  1.72,  2.29,  2.91, 3.65, 4.39, 5.25];   % slip
yITF   = scaling*[0, 10.14, 11.43, 12.29, 12.82, 12.97, 12.44, 11.43, 10.52, 9.42, 8.32, 7.03, 5.84];   % shear stress
yITF_fail = 0;                              % failure shear stress to be imposed for the part that coming out of the matrix
 
% composite and matrix
fc     = 48.8;                              % compressive strength
Ac     = 100*100;                           % trabritary area of composite per one yarn
%rho    = 0.01;                             % reinforcement ratio Af/Ac
Am     = (1-rho)*Ac;                        % trabritary area of matrix per one yarn
Em     = 4735*fc^0.5;                       % young modulus of matrix
 
% longitudinal yarn model (yarn ~ fiber)
nf     = 1;                                 % number of fiber
Af     = rho*Ac;                            % fiber area
df     = (4*Af/pi)^0.5;                     % fiber diameter
psi    = nf*pi*df;                          % perimeter of the fiber
xFIB   = [0, 460/180000, 0.010];            % strain
yFIB   = [0, 460,          470];            % stress
maxSfu = 469;                               % maximum fiber stress, if stress>maxStu, terminate program
effEf  = 1.0;                               % Young modulus efficiency factor
 
% spring model simulating transverse yarn junction
locSPR = [];                                % absolute locations of springs, if no spring, set locSPR=[];
xSPR   = [0,1];                             % slip, if no spring, set xSPR=[0,1];
ySPR   = [0,1];                             % spring force, if no spring, set xSPR=[0,1];
ySPR_fail = 0;                              % failure spring force to be imposed for the spring that coming out of the matrix
 
% matrix strength distribution
%stateNo = 4;                               % state number to generate random number; need for ckPat=3,4,5
%avgST   = 4.4;                             % average matrix strength
%stdST   = 0.02;                            % standard deviation of matrix strength
minST   = 0.5*avgST;                        % trim random generagted strength >= minST
maxST   = 1.5*avgST;                        % trim random generagted strength <= maxST
minSpac = 20;                               % minimum crack spacing
oneCk   = 1;                                % detect one crack at a time
manualPickCrack=0;                          % skip manual selecting crack
 
% numerical analysis parameters
kWt_int  = 0.25;                            % weighting for next increment of interface secant modulus
EfWt_int = 0.25;                            % weighting for next increment of fiber secant modulus
gWt_int  = 0.25;                            % weighting for next increment of spring secant modulus
tor      = [0.0001,0.0005,0.001,0.005];     % multi-level of torlerance for checking material convergence (very small....very large)
accLevel = 1;                               % starting level of torelence, tor(accLevel) will be first used
minIter       = 10;                         % If converge with iterations < minIter, decrease accLevel by one
maxIter       = 1000;                       % If not converge with iterations > maxIter, increase accLevel by one
minCount      = 5;                          % If converge with iterations < minCount, increase weighting factors
maxCount      = 10;                         % If not converge with iterations > maxIter, decrease weighting factors
factCutBack   = 0.95;                       % reduced kWt or EfWt or gWt by a factor of factCutBack
factSpeed     = 1.05;                       % increase kWt or EfWt or gWt by a factor of factSpeed
unstableSlip  = 0.1*L;                      % if slip > unstableSlip, unstable structure quit the program
 
% movie option
showOnOff     = 0;                          % show movie, yes=1, no=0
fps           = 2;                          % fram per second
 
% experimental composite tensile stress strain
hasExpr1      = 0;                          % 0=no experimental data, 1=has
fnameExpr1    = '';                         % file name containing experimental data
startLine1    = 2;                          % starting line (based index=1) for each input file specified in fnameArray 
xValCol1      = 1;                          % column number in the input file that contains x values
yValCol1      = 2;                          % column number in the input file that contains y values
bk0tb1cm2_1   = 2;                          % delimiter for each input file (0=' ', 1='/t', 2=',')
 
% experimental crack spacing strain
hasExpr2      = 0;                          % 0=no experimental data, 1=has
fnameExpr2    = '';                         % file name containing experimental data
startLine2    = 1;                          % starting line (based index=1) for each input file specified in fnameArray 
xValCol2      = 2;                          % column number in the input file that contains x values
yValCol2      = 2;                          % column number in the input file that contains y values
bk0tb1cm2_2   = 0;                          % delimiter for each input file (0=' ', 1='/t', 2=',')
 
% printing setup
colWidth      = 12;                         % width of each column
digit         = 4;                          % decimal used

% ------------------------------------
% Part B: Calculation
% ------------------------------------
% slope of material models
kITF = (yITF(2:end)-yITF(1:end-1))./(xITF(2:end)-xITF(1:end-1));
kFIB = (yFIB(2:end)-yFIB(1:end-1))./(xFIB(2:end)-xFIB(1:end-1));
kSPR = (ySPR(2:end)-ySPR(1:end-1))./(xSPR(2:end)-xSPR(1:end-1));

% select the maximum slope as a reference to calculate the fraction change of stiffness  
% residual = (stiffness_curent-stiffness_previous)/ max_slope v.s. torlerance
% this methodology is bettter than using previous stiffness as a reference
% because previous stiffness can be very low in the late post peak, leading
% to unstable algorithm
kITFmax = max(kITF);
kFIBmax = max(kFIB);
kSPRmax = max(kSPR);

% discretization along the length of fiber
h  = L/(nNode-1);                       % segment length
x  = linspace(0,L,nNode);               % nodal locations
ef = linspace(0,0,nNode);               % nodal fiber strain
k  = linspace(kITF(1),kITF(1),nNode);   % nodal stiffness (shear flow-slip) NOT shear stress-slip
Ef = linspace(kFIB(1),kFIB(1),nNode);   % nodal stiffness (fiber stress-strain)

% generate random matrix strength
randn('state', stateNo);
wkLoc  = linspace(0,L,round(L/minSpac)+1);
wkNode = round(wkLoc/h)+1;
prob   = normrnd(avgST,stdST,length(wkNode),1);
for i=1:length(prob)
    if prob(i) > maxST 
        prob(i) = maxST;
    elseif prob(i) < minST 
        prob(i) = minST;
    end
end
mtxST  = linspace(maxST,maxST,nNode);
mtxST(wkNode) = prob;

% ploting material models
% note that later on in the Part C: show results, 
% the first 3 plots must repeat the first 3 figures below: fiber, interface, spring and matrix envelop models
% let users view material models 1,2,3 and cracking criterion 4 before
% continue the analysis
figure(1)
plot(xFIB,yFIB,'-ro'),title('fiber stress strain'), xlabel('strain'), ylabel('stress'); hold on;

figure(2)
plot(xITF,yITF,'-ro'),title('shear strength - slip'), xlabel('slip'), ylabel('shear stress'); hold on;

figure(3)
plot(xSPR,ySPR,'-ro'),title('spring force - slip'), xlabel('slip'), ylabel('spring force'); hold on;

figure(4)
plot(x,mtxST,'-r.'),title('matrix strength - dist'), xlabel('dist'), ylabel('matrix strength'); hold on;
 
disp('press any key to continue');
%pause

N(1,1:nNode) = 1;                       % initial segment 
N(2,1:nNode) = x;                       % x, distanct
N(3,1:nNode) = kITF(1);                 % k, initial stiffness  (shear flow-slip)
N(4,1:nNode) = kFIB(1);                 % Ef, initial stiffness  (shear flow-slip)
N(5,1:nNode) = -1;                      % g, < 0 for no spring and >=0 for having spring
nodeSPR      = round(locSPR/h);
N(5,nodeSPR) = kSPR(1);                 % assign nodal springs at spring location             
N(6,1:nNode) = mtxST;                   % assign nodal springs at spring location  

% initialize the arrays to record material responses at each increment
SA1 = [];   TA1 = [];       % recorded shear stress - slip response
SN1 = [];   SF1 = [];       % recorded fiber stress - strain response
DS1 = [];   GF1 = [];       % recorded spring force - slip response
ec_cspac    = [];
ec_avgstdCS = [];
prevNseg    = 1;

% initialize load end slip response
D  = [0];                   % incremental end slip
P  = [0];                   % incremental pullout load
CW = [0];                   % incremental crack width
strain = D/L;               % strain of a composite
stress = P/Ac;              % stress of a composite
avg_EM = [0];               % average matrix strain over clear span
avg_SM = [0];               % average matrix stress over clear span
avg_SN = [0];               % average fiber strain over clear span
avg_SF = [0];               % average fiber stress over clear span
avg_FM = [0];               % average matrix force over clear span
avg_FF = [0];               % average fiber force over clear span

count    = 0;
selCount = 0;
adding   = 0;
errCode  = 'none';
for i=2:nInc                % start at i=2, i=1, P=D=0
    % set initial
    ckWidth   = 0;          % a variable to collect crack width of multimple segment at each load step
    endSlip   = 0;          % end slip at both ends   
    S  = [];    T  = [];    SN = [];    SF = [];
    FF = [];    FM = [];    SM = [];    EM = [];
    DS = [];    GF = [];    CS = [];
    
    % apply load step
    P(i) = (i-1)/(nInc-1)*totLoad;

    x    = N(2,:);          % locations
    nSeg = N(1,end);        % current number of segments in a specimen  
    for j = 1:nSeg
        bc = 3;             % bc=3 F(1)=P, F(L)=P   (crack segments)

        % trim data to smaller set of array
        seg_ID  = N(1,:) == j;
        segX    = N(2,seg_ID);
        k       = N(3,seg_ID);
        Ef      = N(4,seg_ID);
        gSeg    = N(5,seg_ID);
        
        segL    = segX(end)-segX(1);
        g_ID    = gSeg >= 0;                            % select only the node that have spring g>=0        
        nodeSPR = find(g_ID==1);
        g       = gSeg(nodeSPR);
        
        % maximum slip at each nodal locations before sliding
        maxSlipITF_Rt = segX(end) - segX + h;           % maximum slip for each node
        maxSlipITF_Lf = segX(1)   - segX - h;
        maxSlipSPR_Rt = segX(end) - segX(g_ID) + h;     % maximum slip for each nodal spring
        maxSlipSPR_Lf = segX(1)   - segX(g_ID) - h;     % maximum slip for each nodal spring
        
        % forming coefficient matrix
        flag     = 1;                                   % check unstability, fiber breakage
        converge = 0;
        iter     = 0;
        kWt      = kWt_int;      
        EfWt     = EfWt_int;        
        gWt      = gWt_int;
        while converge == 0 
            iter = iter + 1;
            s    = Linear_Equation(bc,P(i),h,psi,Af,Ef,k,g,nodeSPR);
            sAbs = abs(s);          % absolute slip (positive number)
            if iter > maxIter
                if accLevel < length(tor)
                    accLevel = accLevel + 1;
                    iter = minCount;
                else
                    flag    = 0;        % quit program
                    comment = 'iter > maxIter, for the lowest acceptable accuracy, quit the program !!!';
                    errCode = 'maxIter';
                    break;              % get off while loop
                end
            end
            if max(sAbs) > unstableSlip
                flag    = 0;        % quit program
                comment = 'unstable slip occurs, quit the program !!!';
                errCode = 'unstableSlip';
                break;              % get off while loop
            end

            % update material secant stiffness
            % interface model
            kPrev = k;
            kNext = Secant_Modulus(xITF,yITF,kITF,sAbs);
            for q=1:length(kNext)
                % if slip at the node i greater than its embeded length maxSlipITF, set the kNext to the failure value
                if s(q) >= 0
                    if s(q) > maxSlipITF_Rt(q)
                        kNext(q) = yITF_fail/s(q);
                    end
                else
                    if s(q) < maxSlipITF_Lf(q)
                        kNext(q) = yITF_fail/abs(s(q));
                    end
                end
            end

            resd_k  = max(abs(kNext-kPrev)/kITFmax);
            k       = (1-kWt)*kPrev + kWt*kNext;
            if resd_k > tor(accLevel) && mod(iter,maxCount)==0
                kWt = kWt*factCutBack;
            end
            
            % fiber model
            ef      = [(s(2)-s(1))/h, (s(3:end)-s(1:end-2))/(2*h), (s(end)-s(end-1))/h];
            EfPrev  = Ef;
            EfNext  = Secant_Modulus(xFIB,yFIB,kFIB,ef);
            resd_Ef = max(abs(EfNext-EfPrev)/kFIBmax);
            Ef      = (1-EfWt)*EfPrev + EfWt*EfNext;
            if resd_Ef > tor(accLevel) && mod(iter,maxCount)==0
                EfWt = EfWt*factCutBack;
            end

            % spring model
            resd_g = 0;
            nSPR   = length(nodeSPR);
            if nSPR > 0
                gPrev = g;
                gNext = Secant_Modulus(xSPR,ySPR,kSPR,sAbs(nodeSPR));
                
                for n=1:nSPR
                    node = nodeSPR(n);
                    % if slip at the spring location i greater than its embeded length maxSlipSPR, set the gNext to the failure value
                    if s(node) >= 0
                        if s(node) > maxSlipSPR_Rt(n)
                            gNext(n) = ySPR_fail/s(node);
                        end
                    else
                        if s(node) < maxSlipSPR_Lf(n)
                            gNext(n) = ySPR_fail/abs(s(node));
                        end
                    end
                end

                resd_g = max(abs(gNext-gPrev)/kSPRmax);
                g = (1-gWt)*gPrev + gWt*gNext;
                if resd_g > tor(accLevel) && mod(iter,maxCount)==0
                    gWt = gWt*factCutBack;
                end
            end

            if (resd_k <= tor(accLevel)) && (resd_Ef <= tor(accLevel)) && (resd_g <= tor(accLevel))
                converge = 1;
                if iter < minIter && accLevel > 1
                    accLevel = accLevel - 1;
                end
                if iter < minCount
                    kWt  = kWt*factSpeed;
                    EfWt = EfWt*factSpeed;
                    gWt  = gWt*factSpeed;
                end
            else
                converge = 0;
            end
            disp(strcat('inc = ',num2str(i),'..seg = ',num2str(j),'..iter = ',...
            num2str(iter),'..resd_k = ',num2str(resd_k),'..resd_Ef = ',...
            num2str(resd_Ef),'..resd_g = ',num2str(resd_g)));
        end

        if flag==1
            % update material stiffness
            N(3,seg_ID)   = k;
            N(4,seg_ID)   = Ef;
            gSeg(nodeSPR) = g;
            N(5,seg_ID)   = gSeg;
            
            % pullout load and end slip
            if nSeg == 1
                endSlip = sAbs(1) + sAbs(end);
            else
                if j==1
                    endSlip = endSlip + sAbs(1);
                    ckWidth  = ckWidth + sAbs(end);
                elseif j==nSeg
                    endSlip = endSlip + sAbs(end);
                    ckWidth  = ckWidth + sAbs(1);
                else
                    ckWidth  = ckWidth + sAbs(1) + sAbs(end);
                end
            end
            CS   = [CS,segL];               % colect crack spacing information
            
            % distributed force
            t  = k.*s;              % nodal shear stress  (sign depending of slip)
            sf = ef.*Ef;            % nodal fiber stress
            ff = sf*Af;             % nodal fiber force
            ds = sAbs(g_ID);        % spring slip
            gf = g.*ds;             % spring force
            fm = P(i)-ff;           % nodal concrete force
            sm = fm/Am;             % nodal matrix stress
            em = sm/Em;             % nodal matrix strain                        
            
            % addd vectors of each segment for a whole specimen length
            S  = [S,s];     T  = [T,t];
            SN = [SN,ef];   SF = [SF,sf];
            DS = [DS,ds];   GF = [GF,gf];
            FF = [FF,ff];   FM = [FM,fm];   
            SM = [SM,sm];   EM = [EM,em];                                  
        else
            break;  % get off the seg loop (j)
        end            
    end
    
    if flag == 1
        % record some data point for tracking material response, also include the maximum value to monitor how far it goes
        nNode     = length(S);
        nRecPt_ID = round(linspace(1,nNode,nRecPt));
        SA1       = [SA1,abs(S(nRecPt_ID)),max(abs(S))];
        TA1       = [TA1,abs(T(nRecPt_ID)),max(abs(T))];
        SN1       = [SN1,SN(nRecPt_ID),max(SN)];
        SF1       = [SF1,SF(nRecPt_ID),max(SF)];
        DS1       = [DS1,DS];
        GF1       = [GF1,GF];
              
        % calculate average stress strain for each constitutes      
        % calculate average stress strain for each constitutes for clear span only
        temp1     = N(2,:);
        index     = temp1 >= endL & temp1 <=(endL+measL);
        avg_EM(i) = mean(EM(index));        avg_SM(i) = mean(SM(index));
        avg_SN(i) = mean(SN(index));        avg_SF(i) = mean(SF(index));
        avg_FM(i) = mean(FM(index));        avg_FF(i) = mean(FF(index));

        % total elogation within specimen does not include endSlip because the measurement taken in side
        % the specimen, not from one end to the other end
        D(i) = ckWidth + avg_EM(i)*L;   

        if nSeg > 1 
            CW(i) = ckWidth/(nSeg-1);   % average crack width
        else
            CW(i) = 0;
        end
        strain(i) = D(i)/L;         
        stress(i) = P(i)/Ac;
        
        if nSeg > prevNseg
            crackDetected    = 1;
        else
            crackDetected    = 0;
        end
        prevNseg = nSeg;
       
        if showOnOff == 1
            figure(17);
            plot(N(2,:),SM,'-b',N(2,:),N(6,:),'r.'), title('matrix stress - x'),
            xlabel('dist'), ylabel('stress');
            AX = axis;  axis([AX(1),AX(2),0,maxST]);
            MVF17(i-1) = getframe;

            figure(18);
            plot(N(2,:),SF,'-k'), title('fiber stress - x'),
            xlabel('dist'), ylabel('stress');
            AX = axis;  axis([AX(1),AX(2),0,totLoad/Af]);
            MVF18(i-1) = getframe;
        end
             
        if crackDetected == 1 || i == nInc
            count        = count+1;
            ckInc(count) = i;
            if crackDetected == 1
                ckStatus(count) = 1;
            else
                ckStatus(count) = 0;
            end
            % record response at cracking or at the last load step
            respCell{count}  = [x;S;T;SN;SF;SM;FF;FM]';
            
            % record same strain for each crack spacing
            ec(1:length(CS)) = strain(i);
            ec_cspac         = [ec_cspac; [ec;CS]'];
            
            % record average and stardard crack spacing for each strain
            ec_avgstdCS      = [ec_avgstdCS; [strain(i),mean(CS),std(CS)]];
            
            numSolid = mod(count,length(solidType))+1;  % add 1 to avoid mod(x,y) = 0  
            numDash  = mod(count,length(solidType))+1;  % add 1 to avoid mod(x,y) = 0  

            % plotting ditribution
            figure(4); plot(x,SM,solidType(numSolid,:)), title('matrix stress - x'), xlabel('dist'), ylabel('m'); hold on;      
            %figure(5); plot(x,S,solidType(numSolid,:)), title('slip - x'), xlabel('dist'), ylabel('slip'); hold on;
            %figure(6); plot(x,T,solidType(numSolid,:)), title('shear stress - x'), xlabel('dist'), ylabel('shear stress'); hold on;
            %figure(7); plot(x,SN,solidType(numSolid,:)), title('fiber strain - x'), xlabel('dist'), ylabel('fiber strain'); hold on;
            %figure(8); plot(x,SF,solidType(numSolid,:)), title('fiber stress - x'), xlabel('dist'), ylabel('fiber stress'); hold on;
            %figure(9); plot(x,FF,dashType(numDash,:), x,FM,solidType(numSolid,:)), title('force - x'), xlabel('dist'), ylabel('force'); legend('fiber','concrete'); hold on;
        end

        % keep the results and save as previous
        i_prev  = i;        x_prev  = x;        S_prev  = S;        T_prev  = T;
        SN_prev = SN;       SF_prev = SF;       SM_prev = SM;       FF_prev = FF;
        FM_prev = FM;       CS_prev = CS;
        
        % check ultimate fiber strength (brittle assumption)
        if max(SF) > maxSfu
            flag = 0;
            comment = 'fiber break quit the program'
            errCode = 'fiberBreak';
            break;  % quit if fiber break
        end
        % check if strain reaches the point where new crack forms
        endL   = 0;
        measL  = L;
        N      = Add_Crack_Node(N,SM,maxST,endL,measL,oneCk);
    else        
        break;      % if flag = 0, get off the increment loop (i)
    end
end

if flag == 0
    disp(comment);
    count = count + 1;
    if crackDetected == 1
        ckStatus(count) = 1;
    else
        ckStatus(count) = 0;
    end
    ckInc(count)    = i_prev;
    respCell{count} = [x_prev; S_prev; T_prev; SN_prev; SF_prev; SM_prev; FF_prev; FM_prev]';
    
    ec(1:length(CS_prev)) = strain(i_prev);
    ec_cspac        = [ec_cspac; [ec;CS_prev]'];
    ec_avgstdCS     = [ec_avgstdCS; [strain(i_prev),mean(CS_prev),std(CS_prev)]];
else
    disp('!!  analysis done');
end

if manualPickCrack==1
    beep;
    count     = 0;
    keepInc   = [];
    numckStep = size(respCell,2);
    for n=1:numckStep
        tt= strcat('crack step : ',num2str(n),'/',num2str(numckStep));
        tempMatrix = respCell{n};
        x  = tempMatrix(:,1);
        SM = tempMatrix(:,6);

        figure(20)
        plot(x,SM,'-b'), title(tt), xlabel('distance'), ylabel('stress')

        yesNo = input('accept this load step : yes=1,no=0....=');
        if yesNo == 1
            count   = count + 1;
            keepInc = [keepInc,n];
            tempckInc(count) = ckInc(n);
            tempCell{count}  = respCell{n};
        end
    end
    respCell = tempCell;
    ckInc    = ckInc(keepInc);
    ckStatus = ckStatus(keepInc);
end

% trim out the crack spacing according to selected crack increment and the end increment
n = 0;
for i=1:length(ckInc)
    strainInc = strain(ckInc(i));
    for j=1:size(ec_cspac,1)
        if ec_cspac(j,1)==strainInc
            n = n+1;
            temp_ec_cspac(n,1:2) = ec_cspac(j,1:2);
        end
    end
    for j=1:size(ec_avgstdCS,1)
        if ec_avgstdCS(j,1)==strainInc
            temp_ec_avgstdCS(i,1:3) = ec_avgstdCS(j,1:3);
        end
    end
end
ec_cspac    = temp_ec_cspac;
ec_avgstdCS = temp_ec_avgstdCS;

% reading experimental data 1, load deformation curve
exprX1=0;    % initial value for x1
exprY1=0;    % initial value for y1
if hasExpr1==1
    if bk0tb1cm2_1==0
        data = dlmread(fnameExpr1,'',startLine1,0);     % use ' ' as a delimiter
    elseif bk0tb1cm2_1==1
        data = dlmread(fnameExpr1,'\t',startLine1,0);   % use '\t' as a delimiter
    elseif bk0tb1cm2_1==2
        data = dlmread(fnameExpr1,',',startLine1,0);    % use ',' as a delimiter
    end
    exprX1 = data(:, xValCol1);
    exprY1 = data(:, yValCol1);
end

% reading experimental data 2, crack spacing vs strain 
exprX2=0;    % initial value for x2
exprY2=0;    % initial value for y2
if hasExpr2==1
    if bk0tb1cm2_2==0
        data = dlmread(fnameExpr2,'',startLine2,0);     % use ' ' as a delimiter
    elseif bk0tb1cm2_2==1
        data = dlmread(fnameExpr2,'\t',startLine2,0);   % use '\t' as a delimiter
    elseif bk0tb1cm2_2==2
        data = dlmread(fnameExpr2,',',startLine2,0);    % use ',' as a delimiter
    end
    exprX2 = data(:, xValCol2);
    exprY2 = data(:, yValCol2);
end

% ------------------------------------
% Part C: Show result
% ------------------------------------
figure(1)
plot(SN1,SF1,'c.'),title('fiber stress - strain'), 
xlabel('strain'), ylabel('stress'); hold on;

figure(2)
plot(SA1,TA1,'c.'),title('shear stress - slip'), 
xlabel('absolute slip'), ylabel('absoulute shear stress'); hold on;

figure(3)
plot(DS1,GF1,'c.'),title('spring force - slip'), 
xlabel('absolute slip'), ylabel('spring force'); hold on;

figure(4)
plot(N(2,:),N(6,:),'-bo'),title('matrix strength - dist'), 
xlabel('dist'), ylabel('matrix strength'); hold off;

% figures 5 to 9 are distribtuions of material response

figure(10)
plot(exprX1*L,exprY1*Ac,'bo', D,P,'-k.', D,avg_FM,'-g.', D,avg_FF,'-r.'),title('load displacement response'), 
xlabel('displacement'), ylabel('load'), grid; legend('experiment','composite','matrix','fiber');

figure(11)
plot(exprX1,exprY1,'bo', strain,stress,'-k.'),title('composite stress strain response'),
xlabel('strain'), ylabel('stress'), grid; legend('experiment','composite');

figure(12)
plot(strain,avg_SM,'-ro'), title('average cracked concrete tensile stress'), 
xlabel('strain'), ylabel('stress');

figure(13)
plot(strain,CW,'-ro'),title('average crack width - strain response'), 
xlabel('strain'), ylabel('average crack width'), grid;

figure(14)
plot(exprX2,exprY2,'r.', ec_cspac(:,1),ec_cspac(:,2),'c.', ec_avgstdCS(:,1),ec_avgstdCS(:,2),'-bo'),
title('mean crack spacing- strain response'), legend('experiment curve','simulation data points','simulation curve')
xlabel('strain'), ylabel('crack spacing'), grid;

figure(15)
plot(ec_avgstdCS(:,1),ec_avgstdCS(:,3),'-bo'),title('standard dev of crack spacing- strain response'), 
xlabel('strain'), ylabel('stdev of crack spacing'), grid;

figure(16)
plot(xFIB,yFIB,'-bx', avg_SN,avg_SF,'-ro'),title('fiber tensile stress'), 
xlabel('strain'), ylabel('stress'), grid; legend('fiber model','average cracked steel');

% Play the movie
if showOnOff == 1
    figure(17)
    pause(4);
    movie(MVF17,1,fps);

    figure(18)
    pause(4);
    movie(MVF18,1,fps);
end

% ----------------------------------------------------------------------
% Part D,       Writing results to the output file
% ---------------------------------------------------------------------- 
disp('!!  writing result.....please wait');

textHeading = {'L';'measL';'endL';'totLoad';'nNode';'nInc';'nRecPt';'yITF_fail';'ySPR_fail';'nf';'df';'Ac';'Am';'Em';'fc';'rho';...
               'avgST';'stdST';'minST';'maxST';'stateNo';'minSpac';'oneCk';'kWt_int';'EfWt_int';'gWt_int';...
               'accLevel';'minIter';'maxIter';'minCount';'maxCount';'factCutBack';'factSpeed';'unstableSlip'};
valHeading  = [L,measL,endL,totLoad,nNode,nInc,nRecPt,yITF_fail,ySPR_fail,nf,df,Ac,Am,Em,fc,rho...
               avgST,stdST,minST,maxST,stateNo,minSpac,oneCk,kWt_int,EfWt_int,gWt_int,...
               accLevel,minIter,maxIter,minCount,maxCount,factCutBack,factSpeed,unstableSlip];

textData    = {'strain';'stress';'CW';'D';'P';'P/Af';'avg_FM';'avg_FF';'avg_SM';'avg_SF';...
               'xITF';'yITF';'SA1';'TA1';'xFIB';'yFIB';'SN1';'SF1';'xSPR';'ySPR';'DS1';'GF1';...
               'straincs';'cspac';'strain';'avgCS';'stdCS'};
data        = {strain,stress,CW,D,P,P/Af,avg_FM,avg_FF,avg_SM,avg_SF,...
               xITF,yITF,SA1,TA1,xFIB,yFIB,SN1,SF1,xSPR,ySPR,DS1,GF1,...
               ec_cspac(:,1),ec_cspac(:,2),ec_avgstdCS(:,1),ec_avgstdCS(:,2),ec_avgstdCS(:,3)};

NJtableResult = StoreTable(data);
NJ            = NJtableResult{1};
tableResult   = NJtableResult{2};

% (1.1) print input and output response
fid1 = fopen(fname,'w');
fprintf (fid1,'program version : %50s\n\n',progVersion);
       
fprintf (fid1,'%12s\n','tor(1..n)');
for i=1:length(tor)
    fprintf (fid1,'%12.4g',tor(i));
end
fprintf (fid1,'\n\n');

fprintf (fid1,'%12s\n','ckInc(1..n)');
for i=1:length(ckInc)
    fprintf (fid1,'%12.4g',ckInc(i));
end
fprintf (fid1,'\n\n');
fprintf (fid1,'error code = %12s\n\n',errCode);

% (1.2) print heading of result file
PrintHeading(fid1,textHeading,valHeading,colWidth,digit);
fprintf (fid1,'**********************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************\n\n');

% (1.3) print table of result file
PrintTable(fid1,textData,tableResult,NJ,colWidth,digit);

% (2.1) Additional responses for plotting
% print selected response at selected increment to another text file
fid2   = fopen(fname_ckResp,'w');
fprintf (fid2,'%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s\n\n','lineNo','x','S','T','SN','SF','SM','FF','FM');
lineNo = 2;
for n=1:size(respCell,2)
    tempMatrix   = respCell{n};
    ckStatusText = num2str(ckStatus(n));
    ckIncText    = num2str(ckInc(n));
    dispText     = num2str(D(ckInc(n)));
    loadText     = num2str(P(ckInc(n)));
    heading      = strcat('crack status = ',ckStatusText,'...Increment =',ckIncText,'...displacement = ',dispText,'...load = ',loadText);
    lineNo       = lineNo+1;
    fprintf (fid2,'%1s\n',heading);

    [nRow, nCol] = size(tempMatrix);
    for i=1:nRow
        lineNo = lineNo+1;        
        fprintf (fid2,'%12d,',lineNo);
        for j=1:nCol
            fprintf (fid2,'%12.4g,', tempMatrix(i,j));
        end
        fprintf (fid2,'\n');
    end
    fprintf (fid2,'\n');
    lineNo = lineNo+1;
end

fclose all;
disp('!!  Done with writing result');