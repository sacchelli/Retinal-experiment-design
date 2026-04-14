function [JML] = JML_stochastic_direct(Ybar,CovY,controls,idx,nn,sig2,theta, Delta, Acell, B, C, D, cutoff, depth)
% function [JML,Ep_norm] = JML_stochastic_direct(Ybar,CovY,controls,idx,nn,sig2,theta, Delta, Acell, B, C, D, cutoff, depth)
% JML = - 2*log likelihood (up to a constant), for the stochastic model

p=length(theta);
%[n,~]=size(B);
[q,~]=size(C);
[qN,~]=size(Ybar);
N=qN/q;
A = 0*Acell{1};
for i = 1:p
    A = A+ theta(i) * Acell{i};
end
Sigma = lyap(A,D*D');
Sigmap = sig2 * eye(q); 

Nexp=sum(nn);           % total number of experiments
%Lnn=length(nn);

Sigmainv=Sigma_inv_theta(theta, Delta, Acell, C, N, Sigma, Sigmap, cutoff, depth);
Eta=response_deterministic_recursive1D(controls,idx,nn,theta, Delta, q, Acell, B, C, N);
JML= ((1/Nexp)*sum( (Ybar-Eta).*(Sigmainv*(Ybar-Eta)),1) * nn'  ...
                            - sum(log(eig(Sigmainv)))+ trace_part_JML(Sigmainv, nn, CovY))/ (N*q);
end