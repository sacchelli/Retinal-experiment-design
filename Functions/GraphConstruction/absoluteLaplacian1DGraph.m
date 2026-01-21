function L = absoluteLaplacian1DGraph(N)
% ABSOLUTELAPLACIAN1DGRAPH Constructs the 1D absolute Laplacian matrix
%
% Generates a sparse N x N matrix representing the 1D Laplacian stencil 
% using absolute values for all coefficients. Unlike the standard 
% combinatorial Laplacian (L = D - A), this returns the "absolute" 
% version (L = D + A) where all non-zero entries are positive.
%
% Structure:
%   - Diagonal entries: 2
%   - Off-diagonal entries: 1
%
% Inputs:
%   N - Number of nodes/points
%
% Output:
%   L - Sparse N x N matrix [1 2 1] stencil
    
    % 1D Laplacian stencil
    e = ones(N,1);
    L = spdiags([e 2*e e], -1:1, N, N);
end