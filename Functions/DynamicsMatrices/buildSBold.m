function SBold = buildSBold(F, Sigma, N)
% BUILDSBOLD Constructs the symmetric block Toeplitz covariance matrix
%
% This function builds the full covariance matrix SBold for a sequence of 
% N states. It populates the matrix by computing the state autocovariance 
% blocks R(k) = F^k * Sigma and arranging them into a symmetric block 
% Toeplitz structure.
%
% Inputs:
%   F     - Discrete-time state transition matrix (n x n)
%   Sigma - Steady-state covariance matrix (n x n)
%   N     - Number of time steps (horizon)
%
% Output:
%   SBold - Symmetric block Toeplitz covariance matrix (N*n x N*n)

n = size(F,1);

dSigmak = Sigma;
col = zeros(N*n,n);
col(1:n,:) = dSigmak/2;
for i=2:N
    dSigmak=F*dSigmak;
    col((i-1)*n+1:i*n,:)=dSigmak;
end

mat = zeros(N*n,N*n);

for i=1:N
    mat((i-1)*n+1:end,(i-1)*n+1:i*n) = col(1:(N-i+1)*n,:);
end

SBold=mat+mat';

end