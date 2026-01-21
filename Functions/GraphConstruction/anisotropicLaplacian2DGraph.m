function [Lx,Ly] = anisotropicLaplacian2DGraph(N)
% ANISOTROPICLAPLACIAN2DGRAPH Decomposes the 2D Laplacian into x and y components
%
% This function constructs two separate sparse N^2 x N^2 matrices representing 
% the 1D absolute Laplacian applied independently along the two grid dimensions.
% Summing these outputs (Lx + Ly) would reconstruct the standard 5-point 
% stencil absolute 2D Laplacian.
%
% Inputs:
%   N  - Number of points per dimension (total nodes = N^2)
%
% Outputs:
%   Lx - Sparse matrix representing the [1 2 1] stencil along the x-axis (rows)
%   Ly - Sparse matrix representing the [1 2 1] stencil along the y-axis (columns)

% 1D Laplacian stencil
e = ones(N,1);
T = spdiags([e 2*e e], -1:1, N, N);
I = speye(N);

% Lx captures connections between adjacent nodes in the same row
Lx = kron(I,T);

% Ly captures connections between adjacent nodes in the same column
Ly = kron(T,I);
end