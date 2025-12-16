function SigmaBold = buildSigmaBoldOld2(C1, C2, F, Sigma, Sigmap, N)

% Careful that this is meant to work with C1 and C2, so if I have only C,
% I need C1=C, C2=C/2 or vice versa.
% 
% Retired 23/09/2025

n = size(F,1);

spdiagF = sparse(kron(eye(N), F));

spdiagC1 = sparse(kron(eye(N), C1));

spdiagC2 = sparse(kron(eye(N), C2));

spdiagSigmap = sparse(kron(eye(N), Sigmap));

Sigmal = zeros(n,n*N);

spSigmal = sparse(n*N,n*N);

SigmaAux = zeros(n,n);
for i = 0:N-1
    SigmaAux = Sigma+F*SigmaAux*F';

    Sigmal(:, i*n+(1:n)) = SigmaAux; 

    spSigmal(i*n+(1:n), i*n+(1:n)) = SigmaAux; 
end


Msub = zeros(N*n,N*n);

MsubAux = Sigmal;

for i=1:N-1
    MsubAux = F * MsubAux(:,1:(N-i)*n);
    for j=0:N-i-1
        Msub(j*n+(1:n),(j+i)*n+(1:n))= MsubAux(:,j*n+(1:n));
    end
end

SigmaBold = spdiagC1*(Msub+spSigmal+Msub')*spdiagC2'+ spdiagC2*(Msub+spSigmal+Msub')*spdiagC1' + spdiagSigmap;
end
