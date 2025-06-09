function MCR = update_MCR(SCR, W)

W = W/sum(W);

MCR = sum(SCR.*W);

end

