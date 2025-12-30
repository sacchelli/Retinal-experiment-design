function plotDiagonalMax(A, options)
% PLOTDIAGONALMAX Plots the maximum absolute value on each diagonal of a matrix
%
% Inputs:
%   A       - Input matrix (can be sparse or dense)
%   options - Optional struct with fields:
%             .useAbs (default: true) - Use absolute values
%             .logScale (default: false) - Use log scale for y-axis
%             .title (default: 'Max value on each diagonal')
%
% Example:
%   plotDiagonalMax(SigmaBoldinv);
%   plotDiagonalMax(SigmaBoldinv, struct('logScale', true));

if nargin < 2
    options = struct();
end

useAbs = getfield_default(options, 'useAbs', true);
logScale = getfield_default(options, 'logScale', false);
titleStr = getfield_default(options, 'title', 'Max value on each diagonal');

n = size(A, 1);
numDiagonals = 2*n - 1;
diagonalIndices = -(n-1):(n-1);
maxValues = zeros(numDiagonals, 1);

for i = 1:numDiagonals
    k = diagonalIndices(i);
    diagVals = diag(A, k);
    
    if isempty(diagVals)
        maxValues(i) = 0;
    elseif useAbs
        maxValues(i) = max(abs(diagVals));
    else
        maxValues(i) = max(diagVals);
    end
end

% Plot
figure;
if logScale
    semilogy(diagonalIndices, maxValues, 'LineWidth', 1.5);
    ylabel('Max |value| (log scale)');
else
    plot(diagonalIndices, maxValues, 'LineWidth', 1.5);
    ylabel('Max |value|');
end

xlabel('Diagonal index (0 = main diagonal)');
title(titleStr);
grid on;

% Add vertical line at main diagonal
hold on;
yLimits = ylim;
plot([0 0], yLimits, 'r--', 'LineWidth', 1);
hold off;

end

% Helper function
function val = getfield_default(s, field, default)
    if isfield(s, field)
        val = s.(field);
    else
        val = default;
    end
end