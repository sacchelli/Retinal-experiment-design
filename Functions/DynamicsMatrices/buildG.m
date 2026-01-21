function G = buildG(A, B, Delta, prec)
% BUILDG Discretizes the input matrix using the Variation of Constants formula
%
% Approximates the integral G = integral_{0}^{Delta} exp(A*(Delta-tau))*B dtau
% using a numerical integration scheme with 'prec' sub-steps. This matrix
% maps the continuous input to the state change over one sample period Delta.
%
% Inputs:
%   A     - Continuous-time system matrix (n x n)
%   B     - Continuous-time input matrix (n x m)
%   Delta - Sampling time interval (scalar)
%   prec  - Precision (number of integration steps)
%
% Output:
%   G     - Discrete-time input matrix

dt = Delta/prec;

As = sparse(A);

Gt = zeros(size(B));

for i=1:prec
    Gt = Gt + dt * ( As * Gt + B);
end

G=Gt;

end