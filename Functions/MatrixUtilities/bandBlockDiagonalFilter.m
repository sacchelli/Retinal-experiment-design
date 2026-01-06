function A_band = bandBlockDiagonalFilter(A, depth, q, m)
% BANDDIAGONALFILTER Extracts band diagonal structure from a matrix
%
% Inputs:
%   A         - Input matrix (can be sparse or dense)
%   bandwidth - Number of diagonals to keep on each side of main diagonal
%               For example, bandwidth=1 keeps main, super, and sub diagonals
%
% Output:
%   A_band    - Matrix with only band diagonal elements (sparse format)
%
% Example:
%   A = randn(100, 100);
%   A_band = bandDiagonalFilter(A, 5);  % Keep only 5 diagonals on each side

% Keep only elements where |i - j| <= bandwidth

N1 = size(A,1)/q;
N2 = size(A,2)/m;


J = zeros(N1*q, N2*m);
for i = 1:N1
    for j = 1:N2
        % Check if block (i,j) is within depth blocks of diagonal
        if abs(i - j) <= depth
            J((i-1)*q + (1:q), (j-1)*m + (1:m)) = ones(q, m);
        end
    end
end

A_band = sparse(J.*A);

end