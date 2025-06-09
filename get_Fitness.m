function [Func_Val, Cons_Viol] = get_Fitness(Pop, func_No)


[Func_Val,g,h] = CEC2017(Pop, func_No);

g(g<0) = 0;
v1 = sum(g,2);

h = abs(h);
h(h<=0.0001) = 0;
v2 = sum(h,2);

Cons_Viol = v1 + v2;

end







