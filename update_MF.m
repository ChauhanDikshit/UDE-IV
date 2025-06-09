function MF = update_MF(SF, W)

W = W/sum(W);
W;
SF_sq = SF.*SF;
SF;
SF_sq;
MF = sum(SF_sq.*W)/sum(SF.*W);

%MF = sum(W'.*SF_sq)/sum(W'.*SF);
MF;

end

