function SBoldInv = buildSBoldInv(F, Sigma, N)
% BUILDSBOLDINV Computes the exact inverse of the block Toeplitz covariance matrix
%
% This function constructs the precision matrix (inverse covariance) SBoldInv.
% For a Gaussian Markov process, the inverse is a block tridiagonal matrix.
% The construction uses the discrete-time Lyapunov-like residual:
% SigmaDelta = Sigma - F * Sigma * F'.
%
% Inputs:
%   F     - Discrete-time state transition matrix (n x n)
%   Sigma - Steady-state covariance matrix (n x n)
%   N     - Number of time steps (horizon)
%
% Output:
%   SBoldInv - Sparse block tridiagonal inverse covariance matrix (N*n x N*n)

n = size(F,1);

SigmaDelta = Sigma - F * Sigma * F';
mat = sparse(N*n,N*n);

SigmaInv = inv(Sigma);
SigmaDeltaInv = SigmaDelta\eye(n);


mat(1:n,1:n) = SigmaInv + F' * SigmaDeltaInv * F;

mat(end-n+1:end,end-n+1:end) = SigmaDeltaInv;

Skk = SigmaDeltaInv + F' * SigmaDeltaInv * F;

Skkm1 = -SigmaDeltaInv*F;
Skkp1 = Skkm1';

mat(n+1:2*n,1:n) = Skkm1;
mat(1:n,n+1:2*n) = Skkp1;

for i=2:N-1
    mat((i-1)*n+1:i*n, (i-1)*n+1:i*n) = Skk;

    mat(i*n+1:(i+1)*n, (i-1)*n+1:i*n) = Skkm1;

    mat((i-1)*n+1:i*n, i*n+1:(i+1)*n) = Skkp1;
end

SBoldInv = mat;

end
