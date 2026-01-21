function SigmaBoldInv = buildSigmaBoldInvCutoff(C, F, Sigma, Sigmap, N,depth)
% BUILDSIGMABOLDINVCUTOFF Computes a banded sparse inverse of the output covariance
%
% This function utilizes the Woodbury Identity and the Jain/Meurant algorithm
% to compute the inverse covariance matrix, but restricts the calculation to
% a specified 'depth' around the diagonal to maintain sparsity.
%
% Logic:
% 1. SigmaBold = CBold * SBold * CBold' + SigmapBold.
% 2. SigmaBoldInv = SigmapBoldInv - SigmapBoldInv * Cbold * T^-1 * Cbold' * SigmapBoldInv
%    where T = SBoldInv + CBold' * SigmapBoldInv * CBold is block tridiagonal:
%
%       |D1  -Ak' 0                 (0)       |
%       |-Ak  D2  -Ak' 0                      |
%   T = | 0  -Ak  D3  -Ak'                    |
%       |   . . . . . . . . . . . . . .       |
%       |       . . . . . . . . . . . . . .   |
%       |                  -Ak   D(N-1)   -Ak'|
%       |     (0)                 -Ak      DN |
%
% 3. T^{-1} is computed using the optimized Jain et al. (2006) algorithm.
% 4. The resulting SigmaBoldInv is truncated beyond 'depth' blocks and cast to sparse.
%
% Inputs:
%   C      - Observation matrix (q x n)
%   F      - State transition matrix (n x n)
%   Sigma  - State covariance matrix (n x n)
%   Sigmap - Measurement noise covariance matrix (q x q)
%   N      - Number of time steps (horizon)
%   depth  - Bandwidth (number of off-diagonal blocks to compute)
%
% Output:
%   SigmaBoldInv - Sparse banded inverse output covariance matrix (N*q x N*q)

n = size(F, 1);
q = size(C, 1);

%% Precompute common terms
SigmaDelta = Sigma - F * Sigma * F';
SigmaInv = Sigma \ eye(n);
SigmaDeltaInv = SigmaDelta \ eye(n);
SigmapInv = Sigmap \ eye(q);
SigmapInvC = Sigmap \ C;
CtSigmapInvC = C' * SigmapInvC;
FtSigmaDeltaInvF = F' * (SigmaDeltaInv * F);

% Build A_blocks for T matrix
Ablks = cell(N, 1);
Ablks{1} = SigmaInv + FtSigmaDeltaInvF + CtSigmapInvC;
Ablks{N} = SigmaDeltaInv + CtSigmapInvC;
for i = 2:N-1
    Ablks{i} = SigmaDeltaInv + FtSigmaDeltaInvF + CtSigmapInvC;
end

Ak = SigmaDeltaInv * F;
Bblk = Ak';
BblkInv = Bblk\eye(n);

%% Jain Algorithm - Compute R and S sequences
R = cell(N-1, 1);
S = cell(N-1, 1);

% Forward pass for R
R{1} = Ablks{1} \ Bblk;
for i = 2:N-1
    R{i} = (Ablks{i} - Bblk' * R{i-1}) \ Bblk;
end

% Backward pass for S
S{N-1} = Bblk / Ablks{N};
for i = N-2:-1:1
    S{i} = Bblk / (Ablks{i+1} - S{i+1} * Bblk');
end

%% Compute D_i blocks (diagonal blocks of T^{-1})
D = cell(N, 1);
D{1} = inv(Ablks{1} - Bblk * S{1}');
for i = 1:N-2
    D{i+1} = BblkInv * S{i} * (eye(n) + Bblk' * (D{i} * S{i}));
end
D{N} = Ablks{N} \ (eye(n) + Bblk' * (D{N-1} * S{N-1}));

%% Compute SigmaBoldInv using Woodbury identity with T^{-1} blocks
% SigmaBoldInv = SigmapBoldInv - SigmapBoldInv * C * T^{-1} * C' * SigmapBoldInv
% We compute blocks of: SigmapInv - SigmapInvC * T^{-1}_{i,j} * C' * SigmapInv

SigmaBoldInv = zeros(N*q, N*q);

% Precompute C' * SigmapInv for efficiency
CtSigmapInv = C' * SigmapInv;

mat = zeros(q*N,q*N);

for i = 1:N
    idxI = (i-1)*q+1:i*q;
    
    % Diagonal block: SigmapInv - SigmapInvC * D{i} * CtSigmapInv
    mat(idxI, idxI) = SigmapInv - SigmapInvC * (D{i} * CtSigmapInv);
    
    % Off-diagonal blocks: accumulate T^{-1} block products
    % T^{-1}_{i,j} = D{i} * S{i} * S{i+1} * ... * S{j-1} for j > i
    SigmapInvCTInvBlk = SigmapInvC * D{i};
    

    jMax = min(N,i+depth);
    
    for j = i+1:jMax
        idxJ = (j-1)*q+1:j*q;
                
        % Accumulate product
        SigmapInvCTInvBlk = SigmapInvCTInvBlk * S{j-1};
        
        % Compute off-diagonal blocks of SigmaBoldInv
        blk = - SigmapInvCTInvBlk * CtSigmapInv;
        mat(idxI, idxJ) = blk;
        mat(idxJ, idxI) = blk';  % Use symmetry
    end
end

SigmaBoldInv = sparse(mat);

end