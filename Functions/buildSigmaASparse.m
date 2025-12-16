function Sigma = buildSigmaASparse(A, sigma, Delta, prec)

dt = Delta/prec;

sigma2 = sigma*sigma';

n = size(A,1);
S = zeros(n);

As = sparse(A);

At = transpose(As);

for i=1:prec
    S = S + dt * ( As * S + S * At + sigma2);
end

Sigma=S;

end