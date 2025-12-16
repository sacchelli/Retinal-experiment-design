function Sigma = buildSigma(A, sigma, Delta, prec)

h=Delta/prec;

sigma2=sigma*sigma';


dS=expm(h*A);

dSk=h/3*sigma2;

sumdSk=dSk;

index = 2+2*mod(1:prec,2);
index(end)=1;

dSt = dS';
for i=1:prec
    dSk=dS*dSk*dSt;
    sumdSk=sumdSk+dSk*index(i);
end

Sigma=sumdSk;

end