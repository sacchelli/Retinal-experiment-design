function G = buildG(A, B, Delta, prec)

h=Delta/prec;

dG=expm(h*A);

dGk=h/3*eye(length(A));

sumdGk=dGk;

index = 2+2*mod(1:prec,2);
index(end)=1;


for i=1:prec
    dGk=dGk*dG;
    sumdGk=sumdGk+dGk*index(i);
end

G=sumdGk*B;

end