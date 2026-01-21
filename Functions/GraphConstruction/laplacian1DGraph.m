function L = laplacian1DGraph(N)
% LAPLACIAN1DGRAPH Constructs the standard 1D discrete Laplacian matrix
%
% Generates a sparse N x N matrix representing the 1D negative Laplacian 
% operator. This is the standard discrete second derivative (finite difference) 
% with a [1 -2 1] stencil.
%
% Structure:
%   - Diagonal: -2
%   - Off-diagonals: 1
%
% Inputs:
%   N - Number of points
%
% Output:
%   L - Sparse N x N tridiagonal Laplacian matrix

% 1D Laplacian stencil
e = ones(N,1);
L = spdiags([e -2*e e], -1:1, N, N);
end