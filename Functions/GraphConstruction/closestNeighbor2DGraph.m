function L = closestNeighbor2DGraph(N)
    % laplacian2DGraph  Constructs the 2D Laplacian matrix (5-point stencil)
    %   N: number of interior points per dimension
    %   Returns L: sparse (N^2 x N^2) Laplacian matrix 
    
    % 1D Laplacian stencil
    e = ones(N,1);
    T = spdiags([e 0*e e], -1:1, N, N);
    I = speye(N);

    % 2D Laplacian = kron(I,T) + kron(T,I)
    L = (kron(I,T) + kron(T,I)) ;
end