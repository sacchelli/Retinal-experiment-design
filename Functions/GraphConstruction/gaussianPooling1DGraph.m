function M = gaussianPooling1DGraph(N,delta,sigma)
    % 
    %   N: number of points per dimension
    %   
    
    
    e = ones(N,1);
    e1 = exp(-delta^2*e/(2*sigma^2));
    M = diag(e);
    for i = 1:N
        ei = e1.^(i^2);
        ei = ei(1:N-i);
        M = M + diag(ei,i) + diag(ei,-i);
    end
end