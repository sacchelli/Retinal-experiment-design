function value = trAB(A, B)
% TRAB computes the trace of the product A*B.
%   value = TRAB(A, B) calculates sum(sum(A .* B.')) using vectorization.
%   
%   We transpose B first to align B's rows with A's columns in memory,
%   then use a dot product.

    % 1. Transpose B so its memory aligns with A.
    %    (Use .' for non-conjugate transpose to handle complex numbers correctly)
    Bt = B.';
    
    % 2. Compute the dot product of the linearized vectors.
    value = Bt(:).' * A(:);
end