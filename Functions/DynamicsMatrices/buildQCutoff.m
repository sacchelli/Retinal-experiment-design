function Q = buildQCutoff(SigmaBoldInv, dSb, N, depth)
% BUILDQCUTOFF Computes the Fisher Information Matrix using a banded block approximation
%
% This function calculates the p x p information matrix Q using a memory-efficient
% tridiagonal super-block approximation. Instead of forming the full Nq x Nq 
% product matrices, it only computes and stores the diagonal, upper, and 
% lower super-blocks to evaluate the trace: 0.5 * trace(Q_i * Q_j).
%
% This "cutoff" approach is designed for systems where the interaction between
% distant time steps is negligible, defined by the 'depth' parameter.
%
% Inputs:
%   SigmaBoldInv - Inverse of the covariance/system matrix (N*q x N*q, sparse)
%   dSb          - Cell array of length p containing derivatives of the 
%                  covariance matrix with respect to each parameter
%   N            - Number of time steps
%   depth        - Block bandwidth (size of one super-block)
%
% Output:
%   Q            - Symmetric p x p information matrix (Fisher Information Matrix)


q = size(SigmaBoldInv,1)/N;
p = length(dSb); 
mq = depth * q;
totalSize = size(SigmaBoldInv, 1); % Safer than q*N
nbrSuperBlocks = ceil(N/depth);

halfQDiag = cell(p, nbrSuperBlocks);
halfQUpperDiag = cell(p, nbrSuperBlocks-1);
halfQLowerDiag = cell(p, nbrSuperBlocks-1);

% --- 1. Compute Block Components of Q_i = SigmaInv * dSb_i ---
for idx = 1:p
    fprintf('Parameter %d/%d', idx, p);
    tStepStart = tic;
    
    for i = 1:nbrSuperBlocks
        % Current Block Rows
        rowStart = (i-1)*mq + 1;
        rowEnd   = min(totalSize, i*mq);
        rows     = rowStart:rowEnd;
        
        % Window of Columns (covers Previous, Current, and Next blocks)
        % This is required to capture the full product bandwidth
        colStart = max(1, (i-2)*mq + 1); 
        colEnd   = min(totalSize, (i+1)*mq);
        cols     = colStart:colEnd;
        
        % A. Diagonal Block (Block i, i)
        halfQDiag{idx,i} = full(SigmaBoldInv(rows, cols)) * full(dSb{idx}(cols, rows));
        
        % B. Off-Diagonal Blocks (Block i, i+1 and i+1, i)
        if i < nbrSuperBlocks
            % Define indices for the NEXT block (i+1)
            nextRowStart = i*mq + 1;
            nextRowEnd   = min(totalSize, (i+1)*mq);
            nextRows     = nextRowStart:nextRowEnd; % Fixed typo here
            
            % Upper Block (Block i, i+1)
            % Result Size: (Size of i) x (Size of i+1)
            halfQUpperDiag{idx, i} = full(SigmaBoldInv(rows, cols)) * full(dSb{idx}(cols, nextRows));
            
            % Lower Block (Block i+1, i)
            % Result Size: (Size of i+1) x (Size of i)
            halfQLowerDiag{idx, i} = full(SigmaBoldInv(nextRows, cols)) * full(dSb{idx}(cols, rows));
        end
    end
    fprintf(' done (%.3fs).\n', toc(tStepStart));
end

fprintf('Building Q ... ');
tStepStart = tic;
% --- 2. Compute Fisher Information Matrix: 0.5 * Tr(Q_i * Q_j) ---
Q = zeros(p, p);
for i = 1:p
    for j = 1:i
        s = 0;
        for k = 1:nbrSuperBlocks
            % 1. Diagonal overlap: Tr( A_kk * B_kk )
            % sum(sum(X .* Y')) is a fast way to compute Trace(X*Y)
            s = s + trAB(halfQDiag{i,k}, halfQDiag{j,k});
            
            if k < nbrSuperBlocks
                % 2. Off-Diagonal overlaps
                % Tr( A_{k,k+1} * B_{k+1,k} ) -> Upper * Lower
                % s = s + sum(sum(halfQUpperDiag{i,k} .* halfQLowerDiag{j,k}.'));
                s = s + trAB(halfQUpperDiag{i,k},halfQLowerDiag{j,k});
                
                % Tr( A_{k+1,k} * B_{k,k+1} ) -> Lower * Upper
                s = s + trAB(halfQLowerDiag{i,k},halfQUpperDiag{j,k});
            end
        end
        Q(i,j) = 0.5 * s;
        Q(j,i) = Q(i,j); 
    end
end
tStepEnd = toc(tStepStart);

fprintf(['done (', num2str(tStepEnd,3) ,'s). \n'])
end