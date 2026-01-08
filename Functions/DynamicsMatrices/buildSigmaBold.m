function SigmaBold = buildSigmaBold(C1, C2, F, Sigma, Sigmap, N)
% Careful that this is meant to work with C1 and C2, so if I have only C,
% I need C1=C, C2=C/2 or vice versa.
%

n = size(F,1);
q = size(C1,1);


C1Fpowers = cell(N, 1);
C2Fpowers = cell(N, 1);

C1Fk = C1;
C2Fk = C2;
C1Fpowers{1} = C1Fk;
C2Fpowers{1} = C2Fk;

for i = 2:N
    C1Fk = C1Fk * F;
    C2Fk = C2Fk * F;
    C1Fpowers{i} = C1Fk;
    C2Fpowers{i} = C2Fk;
end


% After computing C1Fpowers and C2Fpowers, precompute products with Sigma

C1FpowersSigma = cell(N, 1);
C2FpowersSigma = cell(N, 1);

for i = 1:N
    C1FpowersSigma{i} = C1Fpowers{i} * Sigma;
    C2FpowersSigma{i} = C2Fpowers{i} * Sigma;
end


% presum{k,l} stores all q×q matrices for time k,l

presum = cell(N, N);
for k = 1:N    
    for l = 1:k
        presum{k,l} = C1FpowersSigma{k} * C2Fpowers{l}' + ...
                        C2FpowersSigma{k} * C1Fpowers{l}';
    end
end



SigmaAux = zeros(q*N, q*N);


elements = cell(N,N);

for k = 1:N
    elements{k,1} = presum{k,1};
end

for l = 2:N
    for k = l:N
        elements{k,l} = elements{k-1,l-1} + presum{k,l};
    end
end


for l = 1:N
    SigmaAux((l-1)*q+1:l*q,(l-1)*q+1:l*q) = elements{l,l} + Sigmap;
    for k = l+1:N
        SigmaAux((k-1)*q+1:k*q,(l-1)*q+1:l*q) = elements{k,l};
        
        SigmaAux((l-1)*q+1:l*q,(k-1)*q+1:k*q) = (elements{k,l})';
    end
end


SigmaBold = SigmaAux;

end