% give expression and network, get the rank Gene
% ex is the expression data
%net is the matrix of network 
% d is the paramter
function [rank]=RankGene(ex,net,d)
ex = abs(ex);
norm_ex = ex/max(ex);
net_Spar = sparse(net);
degrees = sum(net);
% for specical data we will give the no edge gene to one edge 
index = find(degrees == 0);
degrees(index) = 1;
D1 = sparse(diag(1./degrees));
A = eye(size(net_Spar)) - d*(net_Spar'*D1);
b = (1-d)*norm_ex;
rank = A\b;