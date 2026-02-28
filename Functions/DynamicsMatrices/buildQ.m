function Q = buildQ(SigmaBoldInv, dSb, N)

% BUILDQ Computes the Q-section of the Fisher Information Matrix.
%
% Each entry (i,j) is defined by:
%   Q(i,j) = 0.5 * trace( inv(Sigma)*dSigma_i * inv(Sigma)*dSigma_j )
%           = 0.5 * trace( halfQ{i} * halfQ{j} )
% where halfQ{k} = inv(Sigma) * dSigma_k
%
% Inputs:
%   SigmaBoldInv - Inverse of the covariance matrix (N*q x N*q, sparse)
%   dSb          - Cell array of length p, derivatives of Sigma w.r.t. each parameter
%   N            - Number of time steps  (unused here, kept for API compatibility)
%
% Output:
%   Q            - Symmetric p x p Fisher information matrix

p = size(dSb, 1);


halfQ = cell(p, 1);
for idx = 1:p
    fprintf('Parameter %d/%d', idx, p);
    tStepStart = tic;

    halfQ{idx} = SigmaBoldInv * dSb{idx};   

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