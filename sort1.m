function sorted_Id = sort1(Func_Val, Cons_Val)

[~, ~, unique_Func] = unique(Func_Val);
[~, ~, unique_Cons] = unique(Cons_Val);

temp = (length(Func_Val)+1)*unique_Cons + unique_Func;

[~, sorted_Id] = sort(temp);

end