function bound = UpperBoundFace(M)


% Je fais l'hypothèse que M est de la forme p x p x K

p = size(M,1);
K = size(M,3);

minbound = Inf;

for l = 1:K
    Ml = M(:, :, l);

    maxtrace = -Inf;
    
    for k = 1:K
        Mk = M(:, :, k);
        maxtrace = max(maxtrace, trace(inv(Ml)*Mk));
    end

    minbound = min(minbound, log(det(Ml)) - p + maxtrace);

end

bound = minbound;

end