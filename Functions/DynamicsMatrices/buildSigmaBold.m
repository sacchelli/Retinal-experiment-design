function SigmaBold = buildSigmaBold(C1, C2, F, Sigma, Sigmap, N)
% BUILDSIGMABOLD Constructs the symmetric block Toeplitz output covariance matrix
%
% This function builds the full covariance matrix SigmaBold for a sequence 
% of N outputs. It uses a dual-matrix formulation (C1, C2) to compute the 
% cross-covariance blocks and incorporates observation noise (Sigmap).
%
% Note: To represent the standard output covariance C*Sigma*C', this function 
% expects C1 = C and C2 = C/2 (or vice versa) such that the sum of 
% cross-products (C1*X*C2' + C2*X*C1') reconstructs C*X*C'.
%
% Inputs:
%   C1     - First observation matrix (q x n)
%   C2     - Second observation matrix (q x n), typically C1/2
%   F      - Discrete-time state transition matrix (n x n)
%   Sigma  - State covariance matrix (n x n)
%   Sigmap - Observation noise covariance matrix (q x q)
%   N      - Number of time steps (horizon)
%
% Output:
%   SigmaBold - Symmetric block Toeplitz output covariance matrix (N*q x N*q)


n = size(F,1);
q = size(C1,1);



dSigmak = Sigma;
col = zeros(N*q,q);
col(1:q,:) = (C1*dSigmak*C2'+C2*dSigmak*C1'+Sigmap)/2;
for i=2:N
    dSigmak=F*dSigmak;
    col((i-1)*q+1:i*q,:)=(C1*dSigmak*C2'+C2*dSigmak*C1');
end

mat = zeros(N*q,N*q);

for i=1:N
    mat((i-1)*q+1:end,(i-1)*q+1:i*q) = col(1:(N-i+1)*q,:);
end

SigmaBold=mat+mat';

end