function code = PrintHeading(fid,textData,valData,colWidth,digit)

textStr = strcat('%',num2str(colWidth),'s,');
valStr  = strcat('%',num2str(colWidth),'.',num2str(digit),'g,');

for i=1:size(textData,1)
    fprintf(fid, textStr, textData{i});
end
fprintf(fid, '\n');
for i=1:length(valData)
    fprintf(fid, valStr, valData(i));
end
fprintf(fid, '\n');
