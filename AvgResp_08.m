% Clear all infomation before reading and doing calculation in matlab
clear all;      close all;      clc;

% ----------------------------------------------------------------------
% Part A,       Input data
% ---------------------------------------------------------------------- 
numSample  = 5;                         % number of input files, defined in fnameArray
fnameArray = strvcat('RND_STL_bond13fcr26rho3_Gen01.dat','RND_STL_bond13fcr26rho3_Gen02.dat','RND_STL_bond13fcr26rho3_Gen03.dat','RND_STL_bond13fcr26rho3_Gen04.dat','RND_STL_bond13fcr26rho3_Gen05.dat','RND_STL_bond13fcr26rho3_Gen06.dat','RND_STL_bond13fcr26rho3_Gen07.dat','RND_STL_bond13fcr26rho3_Gen08.dat');
startLine  = [16,16,16,16,16,16,16,16];          % starting line (based index=1) for each input file specified in fnameArray 

xValCol    = [1,1,1,1,1,1,1,1];               % column number in the input file that contains x values
yValCol    = [9,9,9,9,9,9,9,9];               % column number in the input file that contains y values
bktbcm     = [2,2,2,2,2,2,2,2];               % delimiter for each input file (0=' ', 1='/t', 2=',')
resultFname  = 'RND_STL_bond13fcr26rho3_avgSM_strain.dat';
numGenPt     = 500;                     % number of points used to calculate the average curve in interpolation function
minSampleAvg = 1;                       % the miniumum number of samples used in average (defualt = numSample)
startRowResultTable = 16;               % start print result table at line number, !!! make sure the line number starts below the list of input files
sortRawDataNSample  = 0;                % No = 0 (just stack raw data of each sample), Yes=1 (pull data of all sample, then sort by x from min to max)

% paramter to break averaged curve to numGenPt points equally
nEqSpac = 100;                          % number of equal spacing point 
yTox    = 1/1e-3;                       % give some number for the ratio of y to x for one inch on computer screen


% ----------------------------------------------------------------------
% Part B,       Read data
% ----------------------------------------------------------------------
% reading raw data and store it in a table format
disp('Reading input files');
for s=1:numSample
    fname = deblank(fnameArray(s,:));
    disp(['File name : ', fname]);
    startLineBased0 = startLine(s)-1;                   % convert based 1 to based 0 to be used in dlmread function
    if bktbcm(s)==0
        data = dlmread(fname,'',startLineBased0,0);     % use ' ' as a delimiter
    elseif bktbcm(s)==1
        data = dlmread(fname,'\t',startLineBased0,0);   % use '\t' as a delimiter        
    elseif bktbcm(s)==2
        data = dlmread(fname,',',startLineBased0,0);    % use ',' as a delimiter
    end
    tempX = data(:, xValCol(s));
    [xMax(s), indexRawMax(s)] = max(tempX);             % locate xMax and its index
    [xMin(s), indexRawMin(s)] = min(tempX);             % locate xMin and its index    
    xcol = 2*s-1;
    ycol = 2*s;
    for i=1:indexRawMax(s)
        xyRaw(1:indexRawMax(s),xcol) = data(1:indexRawMax(s), xValCol(s));
        xyRaw(1:indexRawMax(s),ycol) = data(1:indexRawMax(s), yValCol(s));
    end
   
    % stack xy data for all samples
    if s==1
        xyRawStack = xyRaw(1:indexRawMax(s), xcol:ycol);
    else
        xyRawStack = [xyRawStack; xyRaw(1:indexRawMax(s), xcol:ycol)];
    end
end

% ----------------------------------------------------------------------
% Part B,       Analysing Data
% ----------------------------------------------------------------------

% Determine the terminated x value in linear interpolation and average response,
% the data beyound this xBreak is ignored.
xMax    = fliplr(sort(xMax));
xBreak  = xMax(minSampleAvg);
xStart  = max(xMin);
xGen    = linspace(xStart, xBreak, numGenPt);
yInterp = zeros(numGenPt, numSample);
disp('Analysing samples');
for s = 1:numSample
    disp(['Sample : ',num2str(s)]);
    xcol   = 2*s-1;
    ycol   = 2*s;
    xySort = sortrows(xyRaw(1:indexRawMax(s),xcol:ycol),1);     % ascending sort data according to x values in column1
    xSort  = xySort(:,1);
    ySort  = xySort(:,2);
    tempX  = xGen(xGen>=min(xSort) & xGen<=max(xSort));
    yInterp(1:length(tempX),s) = linterp(xSort, ySort, tempX)';
end

disp ('!!! Please wait, program is now averaging the responses');
yData         = yInterp(1,:);
tempY         = yData;
nSampleAvg(1) = length(tempY);
yAvg(1)       = mean(tempY);
yMin(1)       = min(tempY);
yMax(1)       = max(tempY);
yStd(1)       = std(tempY);
iLast = length(xGen)
for i=2:iLast
    yData         = yInterp(i,:);
    tempY         = yData(yData~=0);    % select only nonzero y data in the averaging (zero data means the data of some samples doesn't exist)
    nSampleAvg(i) = length(tempY);
    
    if nSampleAvg(i) < minSampleAvg
        iLast = i-1;
        break;
    end
    
    yAvg(i) = mean(tempY);
    yMin(i) = min(tempY);
    yMax(i) = max(tempY);
    yStd(i) = std(tempY);
end
nSampleAvg = nSampleAvg(1:iLast);
xGen = xGen(1:iLast);
yAvg = yAvg(1:iLast);
yMin = yMin(1:iLast);
yMax = yMax(1:iLast);
yStd = yStd(1:iLast);

yAvgMinusyStd   = yAvg - yStd;
yAvgPlusyStd    = yAvg + yStd;
xyInterpAugment = [xGen', yAvg', yMin', yMax', yStd', yAvgMinusyStd', yAvgPlusyStd', nSampleAvg'];

% sort raw data
if sortRawDataNSample  == 1
    xyRawStack = sortrows(xyRawStack,1);    % ascending sort data according to x values in column1
end

% break the xGen-yAvg curve equally in the x-y space
yScale = yAvg/yTox;     % scale magnitude y to magnitude x
ds   = ((xGen(2:end)-xGen(1:end-1)).^2 + (yScale(2:end)-yScale(1:end-1)).^2).^0.5;
s(1) = 0;
for i=1:length(ds)
    s(i+1) = s(i) + ds(i);
end
sGen    = linspace(s(1),s(end),nEqSpac);
xEqSpac = interp1(s,xGen,sGen);
yEqSpac = interp1(s,yAvg,sGen);

% store key analysed data in a table format for printing in output file
N1 = size(xyRawStack,1);
N2 = size(xyInterpAugment,1);
N3 = length(xEqSpac);
NJ = [N1,N1, N2,N2,N2,N2,N2,N2,N2,N2, N3,N3];
for i=1:N1
    Table(i,1:2) = xyRawStack(i,1:2);
end
for i=1:N2;
    Table(i,3:10) = xyInterpAugment(i,1:8);
end
for i=1:N3;
    Table(i,11:12) = [xEqSpac(i),yEqSpac(i)];
end




% ----------------------------------------------------------------------
% Part C,       Ploting Results
% ---------------------------------------------------------------------- 
lineType = ['r','b','m','g','c','y', 'r','b','m','g','c','y']; % preset color

figure(1)
for s=1:numSample
    plot(xyRaw(1:indexRawMax(s),2*s-1), xyRaw(1:indexRawMax(s),2*s), lineType(s)); hold on;
end
Vaxis = axis;   % record the default axis values and applies them to all plots 
grid, axis([0, xBreak, 0, Vaxis(4)]), title('Measured Responses'), xlabel('xRaw'), ylabel('yRaw'), legend(fnameArray);
hold off;

figure(2)
plot(xyRawStack(:,1), xyRawStack(:,2), 'c.', xGen, yAvg, 'r', xGen, yAvgMinusyStd, 'b', xGen, yAvgPlusyStd, 'm'), 
grid, axis([0, xBreak, 0, Vaxis(4)]), title('Averaged Response'), xlabel('x'), ylabel('y'), legend('measured', 'avg', 'avg-1std', 'avg+1std');

figure(3)
plot(xGen, yAvg, '-r.', xEqSpac, yEqSpac, '--bo'), 
grid, axis([0, xBreak, 0, Vaxis(4)]), title('Averaged Response'), xlabel('x'), ylabel('y'), legend('avg', 'equal spacing');

% ----------------------------------------------------------------------
% Part D,       Writing results to the output file
% ---------------------------------------------------------------------- 
fid = fopen(resultFname,'w');

fprintf (fid,'Data Analysis\n');
fprintf (fid,'\n');
fprintf (fid,'number of samples,                  numSample    = %10.4g\n', numSample);
fprintf (fid,'number of points generate,          numGenPt     = %10.4g\n', numGenPt);
fprintf (fid,'mininum samples used in averaging,  minSampleAvg = %10.4g\n', minSampleAvg);
fprintf (fid,'\n');
fprintf (fid,'sample file names\n');
for s=1:numSample
    fprintf (fid,'%5d,%15s\n', s, fnameArray(s,:));
end

if sortRawDataNSample ==1
    fprintf (fid,'sort option: pull raw xy data of all samples together then sort them according to x from min to max\n');
else
    fprintf (fid,'sort option: just stack raw xy data of each sample\n');
end

numInfoLine = 8;        % fixed number
numLineFeed = startRowResultTable - (numInfoLine + numSample + 1);
for i=1:numLineFeed
    fprintf (fid,'\n');
end

% start printing header of the table at line no. #35
fprintf (fid,'%5s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s\n','i', 'x raw', 'y raw', 'xGen', 'yAvg' , 'yMin', 'yMax', 'yStd', 'yAvg-Std', 'yAvg+Std', 'nSampleAvg','xEqSpac','yEqSpac');
for i=1:size(Table,1)
    fprintf (fid,'%5d,',i);
    for j=1:size(Table,2)
        if (i > NJ(j) & Table(i,j) == 0)    % don't print values 0 if they appear at the end of the of the table.
            fprintf (fid,'%12s,', '');
        else
            fprintf (fid,'%12.4g,', Table(i,j));
        end
    end
    fprintf (fid,'\n');
end
fclose(fid);


