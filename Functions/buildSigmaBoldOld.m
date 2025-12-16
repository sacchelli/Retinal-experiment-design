function SigmaBold = buildSigmaBoldTest(C1, C2, F, Sigma, Sigmap, N)

% Careful that this is meant to work with C1 and C2, so if I have only C,
% I need C1=C2=C/2
% 


n = size(F,1);
%q = size(C1,1);

SigmaAux = Sigma;
Sigmal = Sigma;

diagF = F;
diagC1 = C1;
diagC2 = C2;
diagSigmap = Sigmap;



for i = 2:N
    SigmaAux = Sigma+F*SigmaAux*F';
    %display(SigmaAux)
    Sigmal = blkadd(Sigmal, SigmaAux);

    diagF = blkadd(diagF, F);

    diagC1 = blkadd(diagC1, C1);

    diagC2 = blkadd(diagC2, C2);

    diagSigmap = blkadd(diagSigmap, Sigmap);
end

diagF = sparse(diagF);

Msub = zeros(N*n,N*n);
MsubAux = sparse(Sigmal);
%MsubAux = Sigmal;


for i=1:N-1
    MsubAux = diagF * MsubAux;
    MsubAux = sparse([zeros(n,N*n);[MsubAux(1:end-n,1:end-n),zeros((N-1)*n,n)]]);
    %MsubAux = [zeros(n,N*n);[MsubAux(1:end-n,1:end-n),zeros((N-1)*n,n)]];
    Msub = Msub+ MsubAux;
end


SigmaBold = diagC1*(Msub+Sigmal+Msub')*diagC2'+ diagC2*(Msub+Sigmal+Msub')*diagC1' + diagSigmap
end