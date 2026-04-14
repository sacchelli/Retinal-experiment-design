function [M] = MFisher(design_model,controls,theta, sig2,Delta, q, Acell, B, C,D, N,cutoff,depth,regul)
% function [M] = MFisher(design_model,controls,theta, sig2,Delta, q, Acell, B, C,D, N,cutoff,depth,regul)
% computes the Fisher information matrices for parameters theta and the
% input signals in controls by expanding the response
%                                ------------
 
[K,~]=size(controls);
M = cell(K, 1);
p=length(theta);
[n,~]=size(Acell{1});
[~,m]=size(B);

Sigmap = sig2 * eye(q);

A = zeros(n,n);
for i = 1:p
    A = A+ theta(i) * Acell{i};
end

F = buildF(A, Delta);

Sigma = lyap(A,D*D');
    % 
if cutoff
    % Use a cutoff accelerate the computations by making the matrix banded:
    SigmaBoldInv = buildSigmaBoldInvCutoff(C, F, Sigma, Sigmap, N, depth);
else
    SigmaBoldInv = buildSigmaBoldInv(C, F, Sigma, Sigmap, N);
end   
if strcmp(design_model,'deterministic')
    SigmaBoldInv_design=(1/sig2)*eye(N*dn,N*dn);
elseif strcmp(design_model,'stochastic')
    SigmaBoldInv_design=SigmaBoldInv;
end

%% PARTIAL DERIVATIVES PRE-COMPUTATIONS %%
Hb = cell(p, 1);
dSb = cell(p, 1);

for i = 1:p   
    Ab = [A, zeros(n,n); Acell{i}, A];
    
    Bb = [B; zeros(n,m)];
    Cb = [zeros(q,n), C];
    Cbp = [C, zeros(q,n)];
    Db = [D, zeros(n,n); zeros(n,n), zeros(n,n)];

    [Fb,Gb] = buildFG(Ab,Bb, Delta);

    Sigmab = lyap(Ab,Db*Db'); 
    
    if cutoff 
        dSb{i} = buildSigmaBoldCutoff(Cb, Cbp, Fb, Sigmab, zeros(q,q), N, depth);
        Hb{i} = buildHCutoff(Cb, Fb, Gb, N, depth);
    else
        dSb{i} = buildSigmaBold(Cb, Cbp, Fb, Sigmab, zeros(q,q), N);
        Hb{i} = buildH(Cb, Fb, Gb, N);
    end 
end

%% COMPUTATION OF MATRIX Q %%
if strcmp(design_model,'deterministic')
    Q=zeros(p,p);
elseif strcmp(design_model,'stochastic')

    if cutoff
        Q = buildQCutoff(SigmaBoldInv_design, dSb, N, depth);
    else
        Q = buildQ(SigmaBoldInv_design, dSb, N, depth);
    end

end

%% FISHER MATRICES CONSTRUCTION %%

% Please notice that the Fisher matrices are rescaled 
Mi = zeros(p,p);
Hju = zeros(N*q, p);
SBiHju = zeros(N*q, p);

%%
timebar= waitbar(0,'M Fisher computations...');
for i = 1:K
    waitbar(i/K,timebar);
    u = controls{i};
    u = u(:);

    for j = 1:p
        Hju(:,j) = Hb{j} * u;
        SBiHju(:,j) = SigmaBoldInv_design * Hju(:,j);
    end

    for j = 1:p
        for k = 1:j
            uLu_jk = (Hju(:,k)' * SBiHju(:,j) + Hju(:,j)' * SBiHju(:,k))/2;
            Mi(j,k) = uLu_jk + Q(j,k);
            Mi(k,j) = uLu_jk + Q(j,k);
        end
    end  

    M{i} = Mi/(N*q); % Rescale of final Fisher information matrix in order to account for influence of N and q

    %--------------------
    % regularisation of M
    M{i}=M{i}+regul*eye(p);
    %--------------------  
end
close(timebar)
end
