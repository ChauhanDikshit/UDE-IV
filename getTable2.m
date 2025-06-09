str1 = '10k_ResultFunc';
Table1 = [];
Table2 = '';
Table3 = [];

for i = 1 : 7
   str2 = int2str(i);
   result = strcat(str1, str2);
   load (result);
   
   temp_Col1 = [i; Result.Best; Result.Median];
   nums_in_str = strcat(num2str(Result.c(1)), num2str(Result.c(2)), num2str(Result.c(3)));
   temp_Col2 = strcat('[',nums_in_str,']');
   temp_Col3 = [Result.v_Bar; Result.Mean; Result.Worst; Result.STD; Result.SR; Result.vio_Bar];
   
   Table1 = [Table1, temp_Col1];
   Table2 = strcat(Table2,temp_Col2,',');
   Table3 = [Table3, temp_Col3];
end

% Create Table
str = 'Results_10kD.csv';
csvwrite(str, Table1);

Table2(end) = '';
fid=fopen(str, 'a');
% fprintf(fid,Table2);
fprintf(fid, '%s\n', Table2);
fclose(fid);

dlmwrite(str, Table3, 'delimiter', ',', '-append');


