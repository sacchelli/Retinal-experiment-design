function bound = UpperBoundFaceAux(M)

% Je fais l'hypothèse que M est de la forme (M1|M2|...|MK)

p = size(M,1);
K = size(M,2)/p;

minbound = Inf;

for l = 1:K
    Ml = M(:, (l-1)*p+1 : l*p);

    maxtrace = -Inf;
    
    for k = 1:K
        Mk = M(:, (k-1)*p+1 : k*p);
        maxtrace = max(maxtrace, trace(inv(Ml)*Mk));
    end

    minbound = min(minbound, log(det(Ml)) - p + maxtrace);

end

bound = minbound;

end