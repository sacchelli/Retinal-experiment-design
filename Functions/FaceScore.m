function [score, weights] = FaceScore(M)

% Je fais l'hypothèse que M est de la forme p x p x K

p = size(M,1);
K = size(M,3);


w0 = ones(1,K)/K;
Mmu = sum(M .* reshape(w0, 1, 1, K), 3);

Mnew = Mmu;
wnew = w0;


maxIter = 200;
maxBisec = 20;


trajectories=[];


for count = 1:maxIter

    
    Mold=Mnew;
    wold=wnew;
    
    %%%% Ascent direction

    ks = 1;
    scoreTr = -Inf;
    for l = 1:K
        
        scoreTrAux = (trace(inv(Mold)*M(:,:,l))-p)/trace((M(:,:,l)-Mold)'*(M(:,:,l)-Mold));
        
        % This is the normaized gradient, careful because there is a risk
        % of dividing by 0. Below is unnormalized.
        
        % scoreTrAux = trace(inv(Mold)*M(:,:,l));

        if scoreTrAux > scoreTr
            ks = l;
            scoreTr = scoreTrAux;
        end
    end 
    
    Mks = M(:,:,ks);
    wks = zeros(1, K);
    wks(ks)=1;
    
    %%%% Adaptive step 
    
    % The adaptive step is computed via bisection in this method.
    
    a = 0;
    b = 1;
    
    % Assumption is that at a=0, the function is poitive and at b=1 it 
    % is negative. At a it makes sense by construction. However it may 
    % be positive at b in some cases. We don't care.
    
    for i = 1:maxBisec
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

    %improvement = [improvement,log(det(Mnew))/log(det(Mold)) ];
    % trajectories = [trajectories;wnew];
    %trajectories = [trajectories;ks];

end

score = log(det(Mnew));
weights = wnew;

% 
% clf
% %plot(trajectories(:))
% plot3(trajectories(:,1),trajectories(:,2),trajectories(:,3))

end