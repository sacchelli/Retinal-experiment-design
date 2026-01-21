function L = absoluteLaplacian2DGraph(N)
% ABSOLUTELAPLACIAN2DGRAPH Constructs the 2D absolute Laplacian matrix
%
% Generates a sparse N^2 x N^2 matrix representing a 2D Laplacian on 
% a Cartesian grid using a 5-point stencil. This uses absolute notation 
% (L = D + A), where all coefficients are positive.
%
% The 2D Laplacian is constructed via the Kronecker sum of 1D Laplacians:
% L = (I ⊗ T) + (T ⊗ I)
%
% Stencil Structure:
%   - Center (diagonal): 4
%   - Immediate neighbors: 1
%
% Inputs:
%   N - Number of points per dimension (resulting in N^2 total nodes)
%
% Output:
%   L - Sparse N^2 x N^2 absolute Laplacian matrix
    
% 1D Laplacian stencil
e = ones(N,1);
T = spdiags([e 2*e e], -1:1, N, N);
I = speye(N);

% 2D Laplacian = kron(I,T) + kron(T,I)
L = (kron(I,T) + kron(T,I)) ;
end