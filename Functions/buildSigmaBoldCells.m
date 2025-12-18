function SigmaBold = buildSigmaBold(C1, C2, F, Sigma, Sigmap, N)
% Careful that this is meant to work with C1 and C2, so if I have only C,
% I need C1=C, C2=C/2 or vice versa.
%
n = size(F,1);
q = size(C1,1);
spdiagSigmap = sparse(kron(eye(N), Sigmap));

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

% 
% presum{k} stores all q×q matrices for time k (d goes from 0 to k-1)
presum = cell(N, 1);
for k = 1:N
    presum{k} = cell(k, 1);  % For k, we need k entries (d = 0 to k-1)
    for d = 0:k-1
        presum{k}{d+1} = C1Fpowers{k} * Sigma * C2Fpowers{k-d}' + ...
                         C2Fpowers{k} * Sigma * C1Fpowers{k-d}';
    end
end

SigmaAux = zeros(q*N, q*N);

% Off-diagonal blocks
for k = 1:N
    for d = 1:k-1
        % Sum over presum{k-(0:(k-d-1))}{d+1}
        blk = zeros(q, q);
        for offset = 0:(k-d-1)
            blk = blk + presum{k-offset}{d+1};
        end
        SigmaAux((k-1)*q + (1:q), ((k-d)-1)*q + (1:q)) = blk;
    end
end

SigmaAux = SigmaAux + SigmaAux';

% Diagonal blocks
for k = 1:N
    % Sum over presum{k-(0:(k-1))}{1}
    blk = zeros(q, q);
    for offset = 0:(k-1)
        blk = blk + presum{k-offset}{1};
    end
    SigmaAux((k-1)*q + (1:q), (k-1)*q + (1:q)) = blk;
end

SigmaBold = SigmaAux + spdiagSigmap;

end