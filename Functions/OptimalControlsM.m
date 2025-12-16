function [outputs,weights,Mout] = OptimalControlsM(Acol, B, C, theta, sigma, Sigmap, Delta, controls)

% I assume the following shapes: 
% Acol of shape n x n x (p+1) with A_0 = Acol(:,:,p+1), A_i = Acol(:,:,i)
% theta of shape 1 x p
% B of shape n x m
% C of shape q x n
% sigma, Sigma of shape n x n
% Sigmap of shape q x q
% controls of shape m x N x nU


[n,m] = size(B);
p = length(theta);
q = size(C,1);
N = size(controls, 2);
nU = size(controls, 3);

prec = 500;   % parameter for the integration by quadrature, must be even


A=Acol(:,:,p+1);
for i = 1:p
    A=A+theta(i)*Acol(:,:,i);
end




%%%% Step 1: Construction of elements


F = buildF(A, Delta);
%G = buildG(A, B, Delta, prec);
Sigma = buildSigma(A, sigma, Delta, prec);
SigmaBold = buildSigmaBold(C/2, C/2, F, Sigma, Sigmap, N);
SigmaBoldinv = inv(SigmaBold);

%%%% Derivatives 

Bb = [B; zeros(n,m)];
Cb = [zeros(q,n), C];
Cbp = [C, zeros(q,n),];
sigmab = [sigma, zeros(n,n); zeros(n,n), zeros(n,n)];

%

Sigmaboldb = zeros(N*q,N*q,p);
Hb = zeros(N*q,N*m,p);


for i = 1:p
    % Temporary elements
    Abi = [A, zeros(n,n); Acol(:,:,i), A];
    Fbi = buildF(Abi, Delta);
    Gbi = buildG(Abi, Bb, Delta, prec);
    Sigmabi = buildSigma(Abi, sigmab, Delta, prec);
    
    % Saved elements
    
    Hb(:,:,i) = buildH(Cb, Fbi, Gbi, N);
    Sigmaboldb(:,:,i) = buildSigmaBold(Cb,Cbp,Fbi,Sigmabi,zeros(q,q),N);
end

%%%% Qs and Ls

Q = zeros(p,p);
L = zeros(N*m,N*m,p,p);

for i = 1:p
    for j = i:p
        Q(i,j) = 0.5*trace(SigmaBoldinv*Sigmaboldb(:,:,i)*SigmaBoldinv*Sigmaboldb(:,:,j));
        Q(j,i) = 0.5*trace(SigmaBoldinv*Sigmaboldb(:,:,i)*SigmaBoldinv*Sigmaboldb(:,:,j));

        L(:,:,i,j) = 0.5*(Hb(:,:,i)'*SigmaBoldinv*Hb(:,:,j)+Hb(:,:,j)'*SigmaBoldinv*Hb(:,:,i));
        L(:,:,j,i) = 0.5*(Hb(:,:,i)'*SigmaBoldinv*Hb(:,:,j)+Hb(:,:,j)'*SigmaBoldinv*Hb(:,:,i));
    end
end 




%%%% Step 2: Collection of Fisher matrices


M = zeros(p,p,nU);
Mline = zeros(nU,p*(p+1)/2);
for k = 1:nU
    u = controls(:,:,k);
    u = u(:);
    counter=1;
    for i = 1:p
        for j = i:p
            M(i,j,k) = u'*L(:,:,i,j)*u + Q(i,j);
            M(j,i,k) = M(i,j,k);
            Mline(k,counter) = M(i,j,k);
            counter = counter +1;
        end
    end 
end


%%%% Step 3: Convex analysis


% We reduce Mline to its effective rank to avoid near-degeneracy

[U,S,V] = svd(Mline, 'econ');       % Compute SVD

tol = 1e-8;                      % Tolerance for small singular values

rankEff = sum(diag(S) > tol*S(1,1));  % Effective rank

if rankEff<size(Mline,2)
    Mline = Mline * V(:,1:rankEff); % Project onto significant subspace
end

% display(rankEff)


ConvexHullGraph = convhulln(Mline,{'Qx','Q12','Qt'});
% ConvexHullGraph = convhulln(Mline,{'Qx','Q12','Qt'});

vertices = unique(ConvexHullGraph(:))';

%%%% Using a priori estimates on the score

maxScoreVertex = -Inf;
maxScorevertexIndex = 0;

for k = vertices
    
    scorek = logdet(M(:,:,k));

    if scorek >= maxScoreVertex 
        maxScoreVertex = scorek;
        maxScorevertexIndex = k;
    end

    %maxScoreVertex = max(maxScoreVertex,logdet(M(:,:,k)));

end


apriorifaces = [];

for i = 1:size(ConvexHullGraph,1)

    score = UpperBoundFace(M(:,:,ConvexHullGraph(i,:)));
    
    if score >=  maxScoreVertex
        apriorifaces = [apriorifaces,i];
    end
end

ConvHullCap=ConvexHullGraph(apriorifaces,:);



%%%% Computing the score of faces on the cap

scoreMax = -Inf;

faceMax = zeros(1,rankEff);

wMax = 0;

if size(ConvHullCap,1) == 0
    
    outputs = controls(:,:, maxScorevertexIndex);
    weights = 1;
    Mout = M(:,:,maxScorevertexIndex)

else
    
    for i = 1:size(ConvHullCap,1)
    
        [score, w] = FaceScoreDual(M(:,:,ConvHullCap(i,:)));
        
        if score >= scoreMax
            scoreMax = score;
            faceMax = ConvHullCap(i,:);
            wMax = w;
        end
        
       
    end

    outputs = controls(:,:,faceMax);
    weights = wMax;
    Mout = M(:,:,faceMax)
end 


% outputs = controls(:,:,faceMax);
% weights = wMax;


end