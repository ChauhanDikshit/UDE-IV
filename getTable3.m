str1 = '20k_ResultFunc';
Table1 = [];
Table2 = '';
Table3 = [];
JJ=1:28;
for i = JJ
   str2 = int2str(i);
%    str3='_d=10';
   result = strcat(str1, str2);
   load (result);
   
%    temp_Col1 = [i; Result.Best; Result.Median];
    temp_Col1 = [i; Result.Best];
   nums_in_str = strcat(num2str(Result.c(1)), num2str(Result.c(2)), num2str(Result.c(3)));
   temp_Col2 = strcat('[',nums_in_str,']');
   temp_Col3 = [Result.v_Bar; Result.Mean; Result.Worst; Result.STD; Result.SR; Result.vio_Bar];
   %temp_Col3 = [Result.Mean; Result.Worst; Result.STD; Result.SR;];
   Table1 = [Table1, temp_Col1];
   Table2 = strcat(Table2,temp_Col2,',');
   Table3 = [Table3, temp_Col3];
end

% Create Table
str = 'Results_20kD_30D.csv';
csvwrite(str, Table1);

Table2(end) = '';
fid=fopen(str, 'a');
% fprintf(fid,Table2);
fprintf(fid, '%s\n', Table2);
fclose(fid);

dlmwrite(str, Table3, 'delimiter', ',', '-append');


