function [score, weights] = FaceScoreDual(Msimpl)

% I assume that Msimpl is of the form p x p x Ksimpl
% This method aims to use the center of faces as possible directions of
% ascent. First it needs to build them, and at the end to recover the
% weights from them

p = size(Msimpl,1);
Ksimpl = size(Msimpl,3);

K=2*Ksimpl;

M = zeros(p,p,K);
M(:,:,1:Ksimpl) = Msimpl;

for i = 1:Ksimpl
    %indices = setdiff(1:Ksimpl, i);
    M(:,:,Ksimpl+i) = sum(Msimpl(:,:,setdiff(1:Ksimpl, i)),3)/(Ksimpl-1);
end



w0 = ones(1,K)/K;
Mmu = sum(M .* reshape(w0, 1, 1, K), 3);

Mnew = Mmu;
wnew = w0;


maxIter = 50;
maxBisec = 10;


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

end

wfin = wnew(1:Ksimpl);

for i = 1:Ksimpl
    wfin(setdiff(1:Ksimpl, i)) = wfin(setdiff(1:Ksimpl, i)) + wnew(Ksimpl+i)/(Ksimpl-1);
end


score = log(det(Mnew));
weights = wfin;


end