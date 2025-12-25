function L = closestNeighbor1DGraph(N)
    % laplacian1DGraph  Constructs the 1D Laplacian matrix (3-point stencil)
    %   N: number of points per dimension
    %   Returns L: sparse (N x N) Laplacian matrix, escept at the diagonal
    
    % 1D Laplacian stencil
    e = ones(N,1);
    L = spdiags([e 0*e e], -1:1, N, N);
end