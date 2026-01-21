function L = adjacencyLaplacian1DGraph(N)
% ADJACENCYLAPLACIAN1DGRAPH Constructs a 1D Laplacian with uniform weights
%
% Generates a sparse N x N matrix representing a 1D grid connectivity 
% where every element in the 3-point stencil is set to 1. This differs 
% from the standard Laplacian as the diagonal does not sum the degrees 
% of the neighbors, but is instead fixed to 1.
%
% Structure:
%   - Diagonal entries: 1
%   - Off-diagonal entries: 1
%
% Inputs:
%   N - Number of points/nodes in the 1D graph
%
% Output:
%   L - Sparse N x N matrix with [1 1 1] stencil
    
% 1D Laplacian stencil
e = ones(N,1);
L = spdiags([e e e], -1:1, N, N);
end