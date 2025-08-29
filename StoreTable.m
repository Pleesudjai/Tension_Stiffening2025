function NJtableResult = StoreTable(data)

NJ = [];
for n=1:size(data,2)
    var = data{n};
    [s1,s2] = size(var);
    if s1*s2==0     % emty matrix [] s1=1,s2=0
        NJ = [NJ,1];
        m = length(NJ);
        tableResult(1,m) = 0;
    else            % regular matrix
        if s1>s2
            nRow = s1;
            nCol = s2;
        else
            nRow = s2;
            nCol = s1;
            var  = var';
        end
        for j=1:nCol
            NJ = [NJ,nRow];
            m = length(NJ);
            for i=1:nRow
                tableResult(i,m) = var(i,j);
            end
        end
    end
end
NJtableResult = {NJ,tableResult};
