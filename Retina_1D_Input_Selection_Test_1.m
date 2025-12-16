addpath './Functions'
clc

tic


rng(1)

%%%% Data // following the reference.

delta = 0.00389; % cm
sigmaBG = 0.15; % cm
sigmaAG = 0.15; % cm

tauB = 57.5139; % ms
tauA = 86.6298; % ms
tauG = 54.8947; % ms

wm = 0.0128109; % kHz
wp = 0.0188413; % kHz
wBG = 0.0131284; % kHz 
wAG = 0.0137853; % kHz 

%%%% Numerical data

dn = 10; % size of one layer

n = 3 * dn;

m = dn;

p = 7;

q = dn;

N = 100; % Time window


%%%% Model

Gamma = laplacian1DGraph(dn);
GammaBG = gaussianPooling1DGraph(dn,delta,sigmaBG);
GammaAG = gaussianPooling1DGraph(dn,delta,sigmaAG);


id = speye(dn);
zz = sparse(dn,dn);

Acol = zeros(n,n,p+1);

Acol(:,:,1) = [-id , zz, zz ; zz, zz, zz ; zz, zz, zz];
Acol(:,:,2) = [zz , zz, zz ; zz, -id, zz ; zz, zz, zz];
Acol(:,:,3) = [zz , zz, zz ; zz, zz, zz ; zz, zz, -id];

Acol(:,:,4) = [zz , -Gamma, zz ; zz, zz, zz ; zz, zz, zz];
Acol(:,:,5) = [zz , zz, zz ; Gamma, zz, zz ; zz, zz, zz];

Acol(:,:,6) = [zz , zz, zz ; zz, zz, zz ; GammaBG, zz, zz];
Acol(:,:,7) = [zz , zz, zz ; zz, zz, zz ; zz, GammaAG, zz];

%Acol(:,:,8) = [zz , zz, zz ; zz, zz, zz ; zz, zz, zz];

Acol(:,:,8) = 1/tauB*Acol(:,:,1) + 1/tauA*Acol(:,:,2) + 1/tauG*Acol(:,:,3) ...
    + wm*Acol(:,:,4) + wp*Acol(:,:,5) + wBG*Acol(:,:,6) + wAG*Acol(:,:,7);





Tf=10*pi;

prec = 500;   % must be even

Delta = Tf/N;



% sigma = 0.5*eye(n);
% Sigmap = 0.5;

sigma = 0.1*eye(n);
Sigmap = 0.1*eye(q);


% sigma = 0.01*eye(n);
% Sigmap = 10^(-4);

B = sparse([eye(dn,m);zeros(2*dn,m)]);

C = sparse([zeros(q,2*dn),eye(q,dn)]);



%%% Contraction of A.

thetaData = [1/tauB, 1/tauA, 1/tauG, wm, wp, wBG, wAG];


thetaTrue = (rand(1,p)-1/2)/2.*thetaData;

% thetaTrue = [0,0,0,0,0,0,0];

A = Acol(:,:,p+1) + sum(Acol(:,:,1:p) .* reshape(thetaTrue, 1, 1, []), 3);
condition = (max(real(eig(A))) < -.25) || (compteur > stop);
compteur = compteur+1;


display(unique(eig(A)))

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
    disp(i)
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
%%

%%%% Collection of controls


%%% Test family 


sdn=sqrt(dn);

dx = 1/sdn;
dy = 1/sdn;




Kbound=100000;
counter = 0;
controls = zeros(m,N,Kbound);


cMax = 1;
lambdaMax = 2*2*pi;

angleNumber = 30;
velocityNumber = 1;
wavelengthNumber = 2 ;

triples = zeros(3,Kbound);

for la = 0 : angleNumber
    for lv = 1 : velocityNumber
        for lw = 0 : wavelengthNumber
            counter = counter+1;
            thetaU = la/angleNumber*pi;
            c = lv/velocityNumber*cMax;
            lambdaU = lw/wavelengthNumber*lambdaMax;
            for k = 1:N
                compteur = compteur+1;
                temp=zeros(sdn,sdn);
                for i = 1:sdn
                    for j = 1:sdn
                        temp(i,j)=1+cos(lambdaU*i*dx*cos(thetaU)+lambdaU*j*dy*sin(thetaU)+k*Delta*c);
                    end
                end
                controls(:,k,counter) = temp(:);
            end
            triples(:,counter) = [thetaU;c;lambdaU];
        end
    end
end

controls = controls(:,:,1:counter);

triples = triples(:,1:counter);


%%

%%%% Collection of Fisher matrices

K = size(controls,3);

M = zeros(p,p,K);

tic

for i = 1:K
    u = controls(:,:,i);
    u = u(:);
    Mi = zeros(p,p);
    L2u = max(1/Tf*u'*u/N,0.000001);
    u = u/sqrt(L2u);

    for j= 1:p
        for k = 1:j
            val = u'*L(:,:,j,k)*u+Q(j,k);
            Mi(j,k) = val;
            Mi(k,j) = val;
        end
    end

    M(:,:,i) = Mi;
    fprintf(['control ',num2str(i),'/',num2str(K),'\n'])
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

    improvement = (logdet(Mnew)-logdet(Mold))/abs(logdet(Mold));

    display(i)
    display(improvement)

end

%%

fprintf(['Max score of a matrix : ', num2str(maxScore), '\n'])

fprintf(['Final score : ', num2str(logdet(Mnew)), '\n'])

fprintf(['Score of averaged matrix : ', num2str(logdet(Minit)), '\n'])





threshold = 0.001;

idx = find(abs(wnew) > threshold);
%fprintf("weights = " + num2str(wnew(idx(1)))+", " + num2str(wnew(idx(2))))
%display(sum(wnew(idx)))
%display(triples(:,idx))

%%
pause(1)
figure(1)
clf
for l=1:N
    for k=1:length(idx)
        subplot(1,length(idx),k)
        
        imagesc(reshape(controls(:,l,idx(k)), sdn, sdn))
        colormap(gray); 
        axis image; 
        title("weight = " + num2str(wnew(idx(k))) + ...
            ", \theta = " + num2str(triples(1,idx(k))/pi) + ...
            "\pi, c = " + num2str(triples(2,idx(k)))+ ...
            ", \lambda = " + num2str(triples(3,idx(k))/pi)+"\pi")
        drawnow
        
        pause(0.02)
    end
end

