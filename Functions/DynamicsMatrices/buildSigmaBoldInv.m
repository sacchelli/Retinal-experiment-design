function SigmaBoldInv = buildSigmaBoldInv(C, F, Sigma, Sigmap, N)
% BUILDSIGMABOLDINV Inverts the output covariance matrix via Woodbury-Jain decomposition
%
% This function provides an efficient O(N) inversion of the Nq x Nq output 
% covariance matrix. It avoids the cubic complexity of direct inversion by 
% exploiting the Markovian structure of the state-space model.
%
% In this algorithm, we compute the inverse of SigmaBold using a
% combination of Woodbury Identity and Jain et al. 2006 paper "Numerically
% Stable Algorithms for Inversion of Block Tridiagonal and Banded Matrices".
% The point is the following:
% SigmaBold = CBold * SBold * CBold' + SigmapBold.
%
% In particular this allows to compute SigmaBoldInv as
% SigmaBoldInv = SigmapBoldInv - SigmapBoldInv * Cbold * ...
% ...inv(SBoldInv + CBold' * SigmapBoldInv *CBold) * Cbold' * SigmapBoldInv
%
% T = SBoldInv + CBold' * SigmapBoldInv * CBold is block tridiagonal, of the
% form
%       |D1  -Ak' 0                 (0)       |
%       |-Ak  D2  -Ak' 0                      |
%   T = | 0  -Ak  D3  -Ak'                    |
%       |   . . . . . . . . . . . . . .       |
%       |       . . . . . . . . . . . . . .   |
%       |                  -Ak   D(N-1)   -Ak'|
%       |     (0)                 -Ak      DN |
%
% T^{-1} is computed using the optimized Jain et al. 2006 algorithm.
%
% Inputs:
%   C      - Observation matrix (q x n)
%   F      - State transition matrix (n x n)
%   Sigma  - Steady-state covariance matrix (n x n)
%   Sigmap - Measurement noise covariance matrix (q x q)
%   N      - Number of time steps (horizon)
%
% Output:
%   SigmaBoldInv - Full inverse output covariance matrix (N*q x N*q)

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
Bblk = Ak'; %%% To respect Jain et al. notation.
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

% D{1} = (Ablks{1} - Bblk * S{1}')\eye(n); % There seem to be some
% regularity issues around this inversion. So we experiment with
% regularization.

TermToInvert = Ablks{1} - Bblk * S{1}';

% Regularization based on the average scale of the matrix elements
% reg_scale = trace(abs(TermToInvert)) / n; 
% epsilon = reg_scale * 1e-12; % Adjust 1e-12 based on precision needs
epsilon = 0; % no regularization.

D{1} = (TermToInvert + epsilon * eye(n)) \ eye(n);

D{2} = BblkInv * S{1} * (eye(n) + Bblk' * ((TermToInvert ) \ S{1}));

for i = 2:N-2
    %D{i+1} = (A_blocks{i+1} - B_block * S{i+1}') \ (eye(n) + B_block' * (D{i} * S{i}));
    D{i+1} = BblkInv * S{i} * (eye(n) + Bblk' * (D{i} * S{i}));
end
D{N} = Ablks{N} \ (eye(n) + Bblk' * (D{N-1} * S{N-1}));

%% Compute SigmaBoldInv using Woodbury identity with T^{-1} blocks
% SigmaBoldInv = SigmapBoldInv - SigmapBoldInv * C * T^{-1} * C' * SigmapBoldInv
% We compute blocks of: SigmapInv - SigmapInvC * T^{-1}_{i,j} * C' * SigmapInv

SigmaBoldInv = zeros(N*q, N*q);

% Precompute C' * SigmapInv for efficiency
CtSigmapInv = C' * SigmapInv;

for i = 1:N
    idxI = (i-1)*q+1:i*q;
    
    % Diagonal block: SigmapInv - SigmapInvC * D{i} * CtSigmapInv
    SigmaBoldInv(idxI, idxI) = SigmapInv - SigmapInvC * (D{i} * CtSigmapInv);
    
    % Off-diagonal blocks: accumulate T^{-1} block products
    % T^{-1}_{i,j} = D{i} * S{i} * S{i+1} * ... * S{j-1} for j > i
    SigmapInvCTinv_block = SigmapInvC * D{i};
    

     
    for j = i+1:N    
        idxJ = (j-1)*q+1:j*q;
        
        
        % Accumulate product
        SigmapInvCTinv_block = SigmapInvCTinv_block * S{j-1};
        
        % Compute off-diagonal blocks of SigmaBoldInv
        blk = -SigmapInvCTinv_block * CtSigmapInv;
        SigmaBoldInv(idxI, idxJ) = blk;
        SigmaBoldInv(idxJ, idxI) = blk';  % Use symmetry
    end
end

end