function bound = LowerBoundFace(M)

% Je fais l'hypothèse que M est de la forme (M1|M2|...|MK)

p = size(M,1);
K = size(M,2)/p;

maxbound = -Inf;

for l = 1:K
    Ml = M(:, (l-1)*p+1 : l*p);
    maxbound = max(maxbound, log(det(Ml)));
end

bound = maxbound;

end