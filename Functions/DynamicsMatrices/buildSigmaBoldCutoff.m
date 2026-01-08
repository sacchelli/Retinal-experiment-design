function SigmaBold = buildSigmaBoldCutoff(C1, C2, F, Sigma, Sigmap, N, depth)
% BUILDSIGMABOLDCUTOFF Builds SigmaBold with band diagonal approximation
%
% Careful that this is meant to work with C1 and C2, so if I have only C,
% I need C1=C, C2=C/2 or vice versa.
%
% Inputs:
%   C1, C2   - Observation matrices
%   F        - State transition matrix
%   Sigma    - Process noise covariance
%   Sigmap   - Measurement noise covariance
%   N        - Number of time steps
%   depth    - Block bandwidth (only compute blocks within depth blocks of diagonal)
%
% Output:
%   SigmaBold - Block covariance matrix with band diagonal structure

n = size(F,1);
q = size(C1,1);

% Precompute powers of F applied to C1 and C2
C1Fpowers = cell(N, 1);
C2Fpowers = cell(N, 1);
C1Fk = C1;
C2Fk = C2;
C1Fpowers{1} = C1Fk;
C2Fpowers{1} = C2Fk;
for i = 2:N
    C1Fk = C1Fk * F;
    C2Fk = C2Fk * F;
    C1Fpowers{i} = C1Fk;
    C2Fpowers{i} = C2Fk;
end

% Precompute products with Sigma for efficiency
C1FpowersSigma = cell(N, 1);
C2FpowersSigma = cell(N, 1);
for i = 1:N
    C1FpowersSigma{i} = C1Fpowers{i} * Sigma;
    C2FpowersSigma{i} = C2Fpowers{i} * Sigma;
end

% Compute presum{k,l} only for pairs within depth of diagonal
presum = cell(N, N);
for k = 1:N
    lMin = max(1, k - depth);
    for l = lMin:k
        presum{k,l} = C1FpowersSigma{k} * C2Fpowers{l}' + ...
                      C2FpowersSigma{k} * C1Fpowers{l}';
    end
end

% Use FULL matrix (not sparse) - faster for moderate N and depth
SigmaAux = zeros(q*N, q*N);

% Build elements cell array using incremental accumulation
% Only allocate cells we'll actually use
elements = cell(N, N);

% Initialize first column (l=1) for ALL k
% This is needed because diagonal blocks accumulate from all presum{k,1}
for k = 1:N
    elements{k,1} = presum{k,1};
end

% Fill elements incrementally - only within depth band from diagonal
for l = 2:N
    kMin = l;  % k >= l (lower triangular)
    kMax = min(N, l + depth);  % Only within depth
    for k = kMin:kMax
        % Use incremental update: elements{k,l} = elements{k-1,l-1} + presum{k,l}
        elements{k,l} = elements{k-1,l-1} + presum{k,l};
    end
end

% Populate SigmaAux - only fill blocks within depth band
for l = 1:N
    % Diagonal block (k=l) - always computed via first column initialization
    SigmaAux((l-1)*q+1:l*q, (l-1)*q+1:l*q) = elements{l,l} + Sigmap;
    
    % Off-diagonal blocks (k > l), only within depth
    kMax = min(N, l + depth);
    for k = l+1:kMax
        SigmaAux((k-1)*q+1:k*q, (l-1)*q+1:l*q) = elements{k,l};
        SigmaAux((l-1)*q+1:l*q, (k-1)*q+1:k*q) = (elements{k,l})';
    end
end

SigmaBold = sparse(SigmaAux);
end