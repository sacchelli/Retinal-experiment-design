function G = buildGASparse(A, B, Delta, prec)

dt = Delta/prec;


n = size(A,1);
S = zeros(n);

As = sparse(A);

Gt = zeros(size(B));

for i=1:prec
    Gt = Gt + dt * ( As * Gt + B);
end

G=Gt;

end