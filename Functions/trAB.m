function value = trAB(A,B)
% Computes the trace of the product A*B, by computing the sum of the elements 
% of the Kronecker product of A and B'.

value = sum(sum(A.*B'));

end