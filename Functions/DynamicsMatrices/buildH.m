function H = buildH(C, F, G, N)
% BUILDH Constructs the full impulse response (Toeplitz) matrix
%
% This function builds the block lower triangular matrix H, representing 
% the mapping from the input sequence to the output sequence over N steps.
% The matrix is composed of Markov parameters: H_k = C * F^(k-1) * G.
%
% Inputs:
%   C - Observation matrix (q x n)
%   F - Discrete-time state transition matrix (n x n)
%   G - Discrete-time input matrix (n x m)
%   N - Prediction horizon / number of time steps
%
% Output:
%   H - Block Toeplitz matrix (N*q x N*m)

q = size(C,1);
m = size(G,2);


dCHk = C;
col = zeros(N*q,m);
col(1:q,:) = dCHk*G;
for i=2:N
    dCHk=dCHk*F;
    col((i-1)*q+1:i*q,:)=dCHk*G;
end

H = zeros(N*q,N*m);

for i=1:N
    H((i-1)*q+1:end,(i-1)*m+1:i*m) = col(1:(N-i+1)*q,:);
end


end