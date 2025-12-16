function [Lx,Ly] = anisotropicLaplacian2DGraph(N)
    
    % 1D Laplacian stencil
    e = ones(N,1);
    T = spdiags([e 2*e e], -1:1, N, N);
    I = speye(N);

    % 2D Laplacian = kron(I,T) + kron(T,I)
    Lx = kron(I,T);
    Ly = kron(T,I);
end