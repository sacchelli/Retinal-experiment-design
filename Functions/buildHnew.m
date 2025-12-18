function H = buildHnew(C, F, G, N)

n = size(C,2);
q = size(C,1);
m = size(G,2);

% CHANGED: Pre-allocate the full result matrix
H = zeros(N*q, N*m);

% CHANGED: Compute and place directly without concatenation
dHk = eye(size(F));
for k = 1:N
    CFkG = C * dHk * G;
    
    % Fill the block diagonal and lower diagonal
    for i = k:N
        row_idx = (i-1)*q + (1:q);
        col_idx = (k-1)*m + (1:m);
        H(row_idx, col_idx) = CFkG;
    end
    
    if k < N
        dHk = dHk * F;
    end
end

end