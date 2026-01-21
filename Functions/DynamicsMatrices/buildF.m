function F = buildF(A, Delta)
% BUILDF Computes the discrete-time state transition matrix
%
% This function converts a continuous-time system matrix A into its 
% discrete-time equivalent F using the matrix exponential (Zero-Order Hold).
%
% Inputs:
%   A     - Continuous-time system matrix (n x n)
%   Delta - Sampling time interval (scalar)
%
% Output:
%   F     - Discrete-time state transition matrix F = exp(A * Delta)

F= expm(Delta*A);
end