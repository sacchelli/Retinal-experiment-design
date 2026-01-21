function Q = buildQ(SigmaBoldInv, dSb, N, depth)

% BUILDQ Computes the Q-section of the Fisher Information Matrix. 
%
% This function calculates a p x p matrix Q where each entry (i,j) is defined
% by 0.5 * trace(inv(Sigma) * dSigma/dp_i * inv(Sigma) * dSigma/dp_j).
% It uses a "Super-Block" approach to exploit the banded structure of the 
% matrices, significantly accelerating the computation for large N.
%
% Inputs:
%   SigmaBoldInv - Inverse of the covariance/system matrix (N*q x N*q, sparse)
%   dSb          - Cell array of length p containing derivatives of the 
%                  covariance matrix with respect to each parameter
%   N            - Number of time steps
%   depth        - Block bandwidth defining the sparsity pattern
%
% Output:
%   Q            - Symmetric p x p information matrix

q = size(SigmaBoldInv,1)/N;
p = size(dSb,1);

halfQ = cell(p, 1);

mq = depth * q; 
nbrSuperBlocks = ceil(N/depth);

for idx = 1:p
    fprintf('Parameter %d/%d', idx, p);
    tStepStart = tic;
    
    % The result will have a bandwidth of 2 super-blocks
    result = zeros(q*N, q*N);
    
    for i = 1:nbrSuperBlocks
        % Super-block row range
        rowStart = (i-1)*mq + 1;
        rowEnd   = min(i*mq, N*q);
        rows = rowStart:rowEnd;
        
        % Because the band is 'depth' (1 super-block), 
        % SigmaInv only has non-zeros in super-columns [i-1, i, i+1]
        kSuperStart = max(1, i-1);
        kSuperEnd   = min(nbrSuperBlocks, i+1);
        kCols = (kSuperStart-1)*mq + 1 : min(kSuperEnd*mq, N*q);
        
        % Extract the 'Super-Strip' (3 super-blocks wide)
        SigmaBoldInvStrip = full(SigmaBoldInv(rows, kCols));
        
        % For the derivative dSb, the same logic applies: 
        % It only has non-zeros in a 3-super-block neighborhood.
        % The result of SigmaBoldInvStrip * dSb_strip will spread across 5 super-blocks
        jSuperStart = max(1, i-2);
        jSuperEnd   = min(nbrSuperBlocks, i+2);
        jCols = (jSuperStart-1)*mq + 1 : min(jSuperEnd*mq, N*q);
        
        % Extract the relevant part of dSb
        % We need the rows that match Sigma's columns (k_cols)
        % and the columns that cover the resulting 5-super-block span
        DSub = full(dSb{idx}(kCols, jCols));
        
        % One high-speed dense multiply
        result(rows, jCols) = SigmaBoldInvStrip * DSub;
    end

    halfQ{idx} = result;
    fprintf(' done (%.3fs).\n', toc(tStepStart));
end

Q=zeros(p,p);

for i = 1:p
    for j = 1:i
        Q(i,j) = 0.5*trAB(halfQ{i}, halfQ{j});
        Q(j,i) = Q(i,j);
    end
end




end