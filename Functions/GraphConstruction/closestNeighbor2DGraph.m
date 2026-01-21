function L = closestNeighbor2DGraph(N)
% CLOSESTNEIGHBOR2DGRAPH Constructs the adjacency matrix for a 2D grid graph.
%   L = CLOSESTNEIGHBOR2DGRAPH(N) returns a sparse (N^2 x N^2) matrix 
%   representing the 4-connectivity (5-point stencil pattern) of a 
%   uniform 2D grid with N nodes per dimension.
%
%   Input:
%       N - Number of nodes along one dimension (scalar)
%
%   Output:
%       L - Sparse adjacency matrix where L(i,j) = 1 if nodes i and j 
%           are neighbors, and 0 otherwise.

% 1D Adjacency stencil (off-diagonals only)
e = ones(N,1);
T = spdiags([e zeros(N,1) e], -1:1, N, N);
I = speye(N);

% 2D Adjacency = kron(I,T) + kron(T,I)
% This connects each node to its North, South, East, and West neighbors.
L = kron(I,T) + kron(T,I);
end