function L = laplacian2DGraph(N)
% LAPLACIAN2DGRAPH Constructs the standard 2D discrete Laplacian matrix
%
% Generates a sparse N^2 x N^2 matrix representing the 2D negative 
% Laplacian operator on a Cartesian grid using a 5-point stencil. 
% This is the standard form (L = A - D) where the sum of each row is zero 
% (excluding boundary nodes).
%
% The 2D Laplacian is constructed via the Kronecker sum:
% L = (I ⊗ T) + (T ⊗ I)
%
% Stencil Structure:
%   - Center (diagonal): -4
%   - Neighbors (N, S, E, W): 1
%
% Inputs:
%   N - Number of interior points per dimension
%
% Output:
%   L - Sparse N^2 x N^2 2D Laplacian matrix

% 1D Laplacian stencil: [1 -2 1]
e = ones(N,1);
T = spdiags([e -2*e e], -1:1, N, N);
I = speye(N);

% 2D Laplacian assembly
L = (kron(I,T) + kron(T,I));
end