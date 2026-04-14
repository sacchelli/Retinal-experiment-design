function J1 = JEfficientStats(theta, y, SigmaBoldHat, controls, idx, nn, Delta, Acell, B, C, N, D, Sigmap)
% 
q = size(C,1); 
p=length(theta);
A = 0*Acell{1};
for i = 1:p
    A = A+ theta(i) * Acell{i};
end
[F,G] = buildFG(A,B, Delta);
H = buildH(C, F, G, N);

Sigma = lyap(A,D*D');

SigmaBoldInv = buildSigmaBoldInv(C, F, Sigma, Sigmap, N); % Assumption is no cutoff, but the function can be amended here

J1 = -logdet(SigmaBoldInv) + trAB(SigmaBoldInv, SigmaBoldHat);

nnTot=sum(nn);
Lnn = length(nn);
for i_input=1:Lnn
    u = controls{idx(i_input)};
    u = u(:);
    diffOutputs = H*u-y(:,i_input);
    J1=J1+nn(i_input)/nnTot*(diffOutputs'*SigmaBoldInv*diffOutputs);
end

end