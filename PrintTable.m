function PrintTable(fid,textData,tableResult,NJ,colWidth,digit)

textStr = strcat('%',num2str(colWidth),'s,');
valStr  = strcat('%',num2str(colWidth),'.',num2str(digit),'g,');

for i=1:size(textData,1)
    fprintf(fid, textStr, textData{i});
end
fprintf(fid, '\n');

[nRows, nCols] = size(tableResult); 
for i=1:nRows
    for j=1:nCols
        if i>NJ(j)
            fprintf (fid, textStr, '');
        else
            fprintf (fid, valStr, tableResult(i,j));
        end
    end
    fprintf (fid,'\n');
end
