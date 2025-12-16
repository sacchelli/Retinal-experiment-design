function SigmaBold = buildSigmaBoldOld3(C1, C2, F, Sigma, Sigmap, N)

% Careful that this is meant to work with C1 and C2, so if I have only C,
% I need C1=C, C2=C/2 or vice versa.
% 

n = size(F,1);

q = size(C1,1);

%spdiagF = sparse(kron(eye(N), F));

%spdiagC1 = sparse(kron(eye(N), C1));

%spdiagC2 = sparse(kron(eye(N), C2));

spdiagSigmap = sparse(kron(eye(N), Sigmap));

Sigmal = zeros(n,n*N);

spSigmal = sparse(q*N,q*N);

SigmaAux = zeros(n,n);
for i = 0:N-1
    SigmaAux = Sigma+F*SigmaAux*F';

    Sigmal(:, i*n+(1:n)) = SigmaAux; 

    spSigmal(i*q+(1:q), i*q+(1:q)) = C1*SigmaAux*C2'+C2*SigmaAux*C1'; 
end


Msub = zeros(N*q,N*q);

MsubAux = Sigmal;

for i=1:N-1
    MsubAux = F * MsubAux(:,1:(N-i)*n);
    for j=0:N-i-1
        Msub(j*q+(1:q),(j+i)*q+(1:q))= C1*MsubAux(:,j*n+(1:n))*C2'+C2*MsubAux(:,j*n+(1:n))*C1';
    end
end

SigmaBold = Msub + spSigmal + Msub' + spdiagSigmap;

end
