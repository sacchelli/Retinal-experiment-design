function H = buildHCutoff(C, F, G, N, depth)
% BUILDHCUTOFF Builds H matrix with band diagonal approximation
%
% Inputs:
%   C     - Observation matrix
%   F     - State transition matrix
%   G     - Input matrix
%   N     - Number of time steps
%   depth - Block bandwidth (only compute blocks within depth blocks of diagonal)
%
% Output:
%   H - Block lower triangular matrix with band structure

n = size(C,2);
q = size(C,1);
m = size(G,2);

% Precompute C*F^k*G for k=0 to min(N-1, depth)
% We only need up to depth terms since blocks further from diagonal are zero
maxK = min(N-1, depth);
CFkG = cell(maxK + 1, 1);
dHk = eye(size(F));
CFkG{1} = C*dHk*G;
for k = 1:maxK
    dHk = dHk*F;
    CFkG{k+1} = C*dHk*G;
end

% Initialize matrix
H = zeros(N*q, N*m);

% Fill only blocks within depth of diagonal
% For column block j (corresponding to time step j), 
% we fill row blocks i where i >= j and i - j <= depth
for j = 1:N
    iMin = j;  % Lower triangular: i >= j
    iMax = min(N, j + depth);  % Only within depth
    
    for i = iMin:iMax
        % Block (i,j) corresponds to C*F^(i-j)*G
        k = i - j;  % This is the power of F
        H((i-1)*q+1:i*q, (j-1)*m+1:j*m) = CFkG{k+1};
    end
end

H = sparse(H);

end