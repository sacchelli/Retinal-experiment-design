function H = buildH(C, F, G, N)

n = size(C,2);
q = size(C,1);
m = size(G,2);

dHk = eye(size(F));
col = [C*dHk*G];
for i=2:N
    dHk=dHk*F;
    col=[col;C*dHk*G];
end

mat=col;

for i=2:N
    col=[zeros(q,m);col(1:end-q,:)];
    mat=[mat,col];
end

H=mat;

end