clear all;
addpath './Functions'
clc

tic


rng(318)

%%%% Data

dn = 10^2; % must be a square for Laplacian2DGraph to work
% dn = 100;
n = 3 * dn;
m = dn;
p = 4;
q = 1;

N = 100;



[Gammax,Gammay] = anisotropicLaplacian2DGraph(sqrt(dn));

% Gamma = diag(ones(dn-1,1),1)+diag(ones(dn-1,1),-1);

id = speye(dn);
zz = sparse(dn,dn);

Acol = zeros(n,n,p+1);

Acol(:,:,1) = [zz , Gammax, zz ; zz, zz, zz ; zz, zz, zz];

Acol(:,:,2) = [zz , Gammay, zz ; zz, zz, zz ; zz, zz, zz];

Acol(:,:,3) = [zz , zz, zz ; Gammax, zz, zz ; zz, zz, zz];

Acol(:,:,4) = [zz , zz, zz ; Gammay, zz, zz ; zz, zz, zz];

Acol(:,:,5) = [-id , zz, zz ; zz, -id, zz ; id, id, -id];





Tf=2*pi;

prec = 500;   % must be even

Delta = Tf/N;



% sigma = 0.5*eye(n);
% Sigmap = 0.5;

sigma = 0.1*eye(n);
Sigmap = 0.1;


% sigma = 0.01*eye(n);
% Sigmap = 10^(-4);


% B = rand(n,m);
% C = randn(q,n);


% B = [ones(dn,m);zeros(2*dn,m)];
% C = [zeros(q,2*dn),ones(q,dn)];

B = [eye(dn,m);zeros(2*dn,m)];
C = [zeros(q,2*dn),randn(q,dn)];

compteur = 0;
condition = false;
stop = 50;
while condition == false
    %Acol = randn(n,n,p+1);
    thetaTrue = [abs(randn(1,4))]/3;
    A = Acol(:,:,p+1) + sum(Acol(:,:,1:p) .* reshape(thetaTrue, 1, 1, []), 3);
    condition = (max(real(eig(A)))<-.25) || (compteur > stop);
    compteur = compteur+1;
end

display(unique(eig(A)))

display(compteur)

display(max(real(unique(eig(A)))))




%%%% Construction of elements


F = buildF(A, Delta);
G = buildGASparse(A, B, Delta,prec);
Sigma = buildSigmaASparse(A, sigma, Delta, prec);
SigmaBold = buildSigmaBold(C,C/2,F,Sigma,Sigmap,N);
SigmaBoldinv = inv(SigmaBold);



Hb = zeros(N*q,N*m,p);
dSb = zeros(N*q,N*q,p);

for i = 1:p
    display(i)
    Ab = [A, zeros(n,n); Acol(:,:,i), A];
    
    Bb = [B; zeros(n,m)];
    Cb = [zeros(q,n), C];
    Cbp = [C, zeros(q,n),];
    sigmab = [sigma, zeros(n,n); zeros(n,n), zeros(n,n)];

    Fb = buildF(Ab, Delta);
    
    Gb = buildGASparse(Ab, Bb, Delta, prec);
    
    Sigmab = buildSigmaASparse(Ab, sigmab, Delta, prec);
    

    Hb(:,:,i) = buildH(Cb, Fb, Gb, N);
    dSb(:,:,i) = buildSigmaBold(Cb,Cbp,Fb,Sigmab,zeros(q,q),N);
    

end



Q = zeros(p,p);
L = zeros(N*m,N*m,p,p);

for i = 1:p
    for j = 1:p
        Q(i,j) = 0.5*trace(SigmaBoldinv*dSb(:,:,i)*SigmaBoldinv*dSb(:,:,j));
        L(:,:,i,j) = (Hb(:,:,i)'*SigmaBoldinv*Hb(:,:,j) + Hb(:,:,j)'*SigmaBoldinv*Hb(:,:,i))/2;
    end
end

toc


%%%% Collection of controls


%%% Test family 


sdn=sqrt(dn);

dx = 1/sdn;
dy = 1/sdn;

Tf=2*pi;

N=100;
Delta = Tf/N;


lambdaU = 1;
vU = 5;



Kbound=100000;
counter = 0;
controls = zeros(m,N,Kbound);

angleNumber = 20;

for l = 0 : angleNumber-1
    counter = counter+1;
    thetaU = l/angleNumber*pi;
    for k = 1:N
        compteur = compteur+1;
        temp=zeros(sdn,sdn);
        for i = 1:sdn
            for j = 1:sdn
                temp(i,j)=1+cos(2*pi*i*dx*cos(thetaU)+2*pi*j*dy*sin(thetaU)+k*Delta*vU);
            end
        end
        controls(:,k,counter) = temp(:);
    end
end

controls = controls(:,:,1:counter);


%%

%%%% Collection of Fisher matrices

K = size(controls,3);

M = zeros(p,p,K);

tic

for i = 1:K
    u = controls(:,:,i);
    u = u(:);
    Mi = zeros(p,p);

    for j= 1:p
        for k = 1:j
            val = u'*L(:,:,j,k)*u+Q(j,k);
            Mi(j,k) = val;
            Mi(k,j) = val;
        end
    end

    M(:,:,i) = Mi;
    fprintf(['control ',num2str(i),'\n'])
end

toc
%%


%%%% Finding the best vertex for initialization

maxScore = -Inf;
maxScoreIndex = 0;

for k = 1:K
    if maxScore < logdet(M(:,:,k))
        maxScore = logdet(M(:,:,k));
        maxScoreIndex = k;
    end
    %maxscorevertex = max(maxscorevertex,log(det(M(:,:,k))));

end

Mtop = M(:,:,maxScoreIndex);

%w = zeros(K,1);
%w(maxScoreIndex) = 1; % weight initialization


maxBisec = 20;

wold = ones(K,1)/K;
wnew = wold;

Minit = sum(M .* reshape(wnew, 1, 1, K), 3);

Mold = Minit;
Mnew = Mold;

%%

maxStepSearch = 100;
i = 0;
improvementThreshold = 0.000001;
improvement = 1;

while (i < maxStepSearch)&& (improvement > improvementThreshold)
    
    i = i+1;

    Mold = Mnew;
    wold = wnew;

    % Step 1: search for highest derivative

    maxDer = -Inf;
    maxDerIndex = 0;
    for k = 1:K                     
        der = trace(inv(Mold)*M(:,:,k));
        if maxDer < der
            maxDer = der;
            maxDerIndex = k;
        end
    end
    
    ks = maxDerIndex;
    
    % Step 2: search for best step

    Mks = M(:,:,ks);
    wks = zeros(K, 1);
    wks(ks)=1;
    
    %%%% Adaptive step 
    
    % The adaptive step is computed via bisection in this method.
    
    a = 0;
    b = 1;
    
    % Assumption is that at a = 0, the function is positive and at b=1 it 
    % is negative. At a it makes sense by construction. However it may 
    % be positive at b in some cases. We don't care.
    
    for j = 1:maxBisec
        c = (a+b)/2;
        if trace(inv((1-c)*Mold+c* Mks)*(Mks-Mold))>=0
            a = c;
        else 
            b = c;
        end
    end
    alpha = (a+b)/2; % Final step
    
    
    Mnew = (1-alpha) * Mold + alpha * Mks;
    wnew = (1-alpha) * wold + alpha * wks;
    
    %display(logdet(Mnew))

    improvement = (logdet(Mnew)-logdet(Mold))/abs(logdet(Mold))

    display(i)
    display(improvement)

end

display(maxScore)

display(logdet(Mnew))

display(logdet(Minit))





threshold = 0.001;

idx = find(abs(wnew) > threshold);
display(wnew(idx))
display(sum(wnew(idx)))

%%
figure(1)
clf
for l=1:N
    for k=1:length(idx)
        subplot(1,length(idx),k)
        imagesc(reshape(controls(:,l,idx(k)), sdn, sdn))
        drawnow
        pause(0.01)
    end
end
