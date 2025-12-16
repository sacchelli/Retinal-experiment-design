function SigmaBold = buildSigmaBoldSparse(C1, C2, F, Sigma, Sigmap, N)

% Careful that this is meant to work with C1 and C2, so if I have only C,
% I need C1=C2=C/2
% 


n = size(F,1);
%q = size(C1,1);

%SigmaAux = Sigma;
SigmaAux = zeros(n,n);

%spdiagF = sparse(n*N,n*N);
spdiagF = sparse(kron(eye(N), F));

%spdiagC1 = sparse(q*N,n*N);
spdiagC1 = sparse(kron(eye(N), C1));

%spdiagC2 = sparse(q*N,n*N);
spdiagC2 = sparse(kron(eye(N), C2));

%spdiagSigmap = sparse(n*q,n*q);
spdiagSigmap = sparse(kron(eye(N), Sigmap));

Sigmal = zeros(n*N,n*N);

%spSigmal = sparse(n*N,n*N);

for i = 0:N-1
    SigmaAux = Sigma+F*SigmaAux*F';
    
    Sigmal(i*n+(1:n), i*n+(1:n)) = SigmaAux; 
end

spSigmal = sparse(Sigmal);

Msub = zeros(N*n,N*n);
%Msub = sparse(N*n,N*n);
MsubAux = spSigmal;



for i=1:N-1
    MsubAux = spdiagF * MsubAux;
    MsubAux = sparse([zeros(n,N*n);[MsubAux(1:end-n,1:end-n),zeros((N-1)*n,n)]]);
    Msub = Msub + MsubAux;
end


SigmaBold = spdiagC1*(Msub+spSigmal+Msub')*spdiagC2'+ spdiagC2*(Msub+spSigmal+Msub')*spdiagC1' + spdiagSigmap;
end