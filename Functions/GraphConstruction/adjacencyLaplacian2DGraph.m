function L = adjacencyLaplacian2DGraph(N)
% ADJACENCYLAPLACIAN2DGRAPH Constructs a 2D Laplacian with unit diagonal
%
% Generates a sparse N^2 x N^2 matrix representing a 2D grid where the 
% self-loop (diagonal) and the neighborhood connections (off-diagonals) 
% are balanced. 
%
% This is constructed via the Kronecker sum of 1D [1, 0.5, 1] stencils, 
% resulting in a 5-point 2D stencil:
%   - Center (diagonal): 1.0 (from 0.5 + 0.5)
%   - Neighbors (N, S, E, W): 1.0
%
% Inputs:
%   N - Number of points per dimension
%
% Output:
%   L - Sparse N^2 x N^2 matrix with unit-weight 5-point stencil

% 1D Laplacian stencil (Diagonal = 0.5 to ensure 2D Diagonal = 1.0)
e = ones(N,1);
T = spdiags([e e/2 e], -1:1, N, N);
I = speye(N);

% 2D Laplacian = kron(I,T) + kron(T,I)
L = (kron(I,T) + kron(T,I));
end