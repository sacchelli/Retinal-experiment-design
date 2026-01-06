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
%   SigmaBold - Block covariance matrix with band diagonal structure (sparse)

n = size(F,1);
q = size(C1,1);

% Measurement noise contribution
spdiagSigmap = sparse(kron(eye(N), Sigmap));

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

% Compute presum terms, but ONLY for d <= depth
% presum{k}{d+1} is only computed if d <= depth
presum = cell(N, 1);

for k = 1:N
    % Only compute up to min(k-1, depth) terms
    maxD = min(k-1, depth);
    presum{k} = cell(maxD + 1, 1);  % d goes from 0 to maxD
    
    for d = 0:maxD
        presum{k}{d+1} = C1FpowersSigma{k} * C2Fpowers{k-d}' + ...
                         C2FpowersSigma{k} * C1Fpowers{k-d}';
    end
end

% Initialize sparse matrix
SigmaAux = sparse(q*N, q*N);

% Off-diagonal blocks: only fill blocks within depth blocks of diagonal
for k = 1:N
    maxD = min(k-1, depth);  % Only process d <= depth
    
    for d = 1:maxD
        % For block (k, k-d), we need to sum presum{k-offset}{d+1}
        % But we only have presum terms up to depth
        % The offset goes from 0 to min(k-d-1, depth-d)
        
        blk = zeros(q, q);
        maxOffset = k - d - 1;  % Maximum possible offset
        
        for offset = 0:maxOffset
            kIdx = k - offset;
            % Only accumulate if this presum term exists (d <= depth for that k)
            if d <= min(kIdx-1, depth) && kIdx <= N
                blk = blk + presum{kIdx}{d+1};
            end
        end
        
        SigmaAux((k-1)*q + (1:q), ((k-d)-1)*q + (1:q)) = blk;
    end
end

% Make symmetric
SigmaAux = SigmaAux + SigmaAux';

% Diagonal blocks
for k = 1:N
    blk = zeros(q, q);
    
    % Sum over presum{k-offset}{1} where offset goes from 0 to k-1
    % But we only have presum{j}{1} for j where 0 <= depth (always satisfied for d=0)
    for offset = 0:(k-1)
        kIdx = k - offset;
        if kIdx >= 1 && kIdx <= N
            blk = blk + presum{kIdx}{1};
        end
    end
    
    SigmaAux((k-1)*q + (1:q), (k-1)*q + (1:q)) = blk;
end

% Add measurement noise
SigmaBold = SigmaAux + spdiagSigmap;

end