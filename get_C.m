function c = get_C(vec, func_No)

[~, g, h] = CEC2017(vec, func_No);

cons_viol = [g,h];
c = zeros(1,3);

for i = 1:size(cons_viol,2)
	if cons_viol(i)>1
		c(1) = c(1) +1;
	elseif cons_viol(i)>0.01
		c(2) = c(2) +1;
	elseif cons_viol(i)>0.0001
		c(3) = c(3) +1;
	end
end

end

