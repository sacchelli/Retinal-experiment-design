function SigmaBold = buildSigmaBoldCutoff(C1, C2, F, Sigma, Sigmap, N, depth)
% BUILDSIGMABOLDCUTOFF Builds a banded approximation of the output covariance matrix
%
% Constructs a symmetric block-banded Toeplitz matrix SigmaBold. Blocks further 
% than 'depth' steps from the diagonal are assumed to be zero, significantly 
% reducing memory usage for large N.
%
% Note: To reconstruct C*Sigma*C', use C1 = C and C2 = C/2. The function 
% calculates (C1*X*C2' + C2*X*C1') to ensure symmetry during construction.
%
% Inputs:
%   C1        - First observation matrix (q x n)
%   C2        - Second observation matrix (q x n), typically C1/2
%   F         - Discrete-time state transition matrix (n x n)
%   Sigma     - Steady-state covariance matrix (n x n)
%   Sigmap    - Observation noise covariance (q x q)
%   N         - Number of time steps (horizon)
%   depth     - Block bandwidth; only terms up to lag 'depth' are computed
%
% Output:
%   SigmaBold - Sparse symmetric block-banded covariance matrix (N*q x N*q)

q = size(C1,1);

maxK = min(N-1, depth);

dSigmak = Sigma;
col = zeros(N*q,q);
col(1:q,:) = (C1*dSigmak*C2'+C2*dSigmak*C1'+Sigmap)/2;
for i=2:maxK
    dSigmak=F*dSigmak;
    col((i-1)*q+1:i*q,:)=(C1*dSigmak*C2'+C2*dSigmak*C1');
end

mat = zeros(N*q,N*q);

for i=1:N
    mat((i-1)*q+1:end,(i-1)*q+1:i*q) = col(1:(N-i+1)*q,:);
end

SigmaBold=sparse(mat+mat');

end




