function A_band = bandDiagonalFilter(A, bandwidth)
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
A_band = triu(tril(A, bandwidth), -bandwidth);

% Ensure output is sparse if input was sparse
if issparse(A)
    A_band = sparse(A_band);
end

end