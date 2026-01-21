function L = closestNeighbor1DGraph(N)
% CLOSESTNEIGHBOR1DGRAPH Constructs the 1D adjacency matrix (off-diagonal only)
%
% This function generates a sparse N x N matrix representing the connectivity 
% of a 1D line graph. It uses a [1 0 1] stencil, meaning it captures 
% connections to the immediate left and right neighbors while leaving the 
% diagonal as zero.
%
% Inputs:
%   N - Number of nodes in the graph
%
% Output:
%   L - Sparse N x N adjacency matrix

% 1D stencil: [1 0 1]
e = ones(N,1);
L = spdiags([e 0*e e], -1:1, N, N);
end