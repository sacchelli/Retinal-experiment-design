function [F,G] = buildFG(A, B, Delta)
% BUILDFG computes the transition matrix and the input matrix using 
% the Variation of Constants formula.
%
% In order to compute F = exp(Delta*A) and the integral 
% G = integral_{0}^{Delta} exp(A*(Delta-t))*B dt using a matrix
% exponentiation trick centered around the matrix exponential of a welll
% chosen block matrix.
%
% Inputs:
%   A     - Continuous-time system matrix (n x n)
%   B     - Continuous-time input matrix (n x m)
%   Delta - Sampling time interval (scalar)
%
% Output:
%   F     - Transition matrix
%   G     - Discrete-time input matrix

[n,m] = size(B);

M = [Delta*A, Delta*B; 
         zeros(m, n+m)];

expM = expm(M);

F=expM(1:n,1:n);
G=expM(1:n,n+1:n+m);

end