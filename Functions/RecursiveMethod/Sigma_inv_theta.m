function SigmaBoldInv = Sigma_inv_theta(theta, Delta, Acell, C, N, Sigma, Sigmap, cutoff, depth)
% function SigmaBoldInv = Sigma_inv_theta(theta, Delta, Acell, C, N, Sigma, Sigmap, cutoff, depth)
% 

p=length(theta);
A = 0*Acell{1};
for i = 1:p
    A = A+ theta(i) * Acell{i};
end
F = buildF(A, Delta);
if cutoff
    % Use a cutoff accelerate the computations by making the matrix banded:
    SigmaBoldInv = buildSigmaBoldInvCutoff(C, F, Sigma, Sigmap, N, depth);
else
    SigmaBoldInv = buildSigmaBoldInv(C, F, Sigma, Sigmap, N);
end  

end