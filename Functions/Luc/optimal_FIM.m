function [Mopt,wopt,indices] = optimal_FIM(MM,winit,crit,KK,maxBisec,eff,remove_inessential,period)
% function [Mopt,wopt,indices] = optimal_FIM(MM,winit,crit,KK,maxBisec,eff,remove_inessential,period)
% Computes the optimal information matrix for the criterion crit in the
% form of a weighted sum of p*p matrices in MM (there are K such matrices)
% MM = cell(K, 1)
% winit = initial vector of weights of length K 
%           (if empty, each matrix in MM receives the weight 1/K)
% crit = Dopt or Lopt
%   For crit = Dopt, log(det(M)) is maximised
%   For crit = Lopt, trace(KK'*M^(-1)*KK) is minimised, with KK a p*m matrix 
%       (possibly a column vector of length p); KK may be empty if crit = Dopt
%   For crit = Dopt and remove_inessential is true, inessential points are
%       removed using the sieve of [Harman & Trnovska, Mathematica Slovaka, 2009]
% The efficiency reached by the algorithm is shown every period iterations 
%
% indices refer to the matrices in MM, wopt is the vector of optimal
% weights (>=0 and summing to 1) and has the same length as indices
% (possibly less than K if crit = Dopt and remove_inessential is true)
% The optimisation uses Fedorov's algorithm with optimal step-size found by
% bisection (using maxBisec iterations -- 20 is largely enough)
% eff in (0,1) = efficiency required with respect to the unknown optimal matrix M*
%   constructed with MM
%   -- if crit = Dopt, (det(Mopt)/det(M*))^(1/p) > eff 
%   -- if crit = Lopt, trace(KK'*M*^(-1)*KK)/trace(KK'*Mopt^(-1)*KK) > eff
% The algorithm may have troubles if crit = Lopt, rank(KK)<p and some
% matrices in MM don't have full rank

[K,~]=size(MM);
[p,~]=size(MM{1});
if isempty(winit)
    winit = ones(K,1)/K;
end
wold = winit;
wnew = wold;
% Compute weighted sum
    Minit = zeros(p,p);
    for k = 1:K
        Minit = Minit + wnew(k) * MM{k};
    end
Mold = Minit;
Mnew = Mold;
indices=1:K; % indices of useful inputs (used in case some matrices are removed from MM 
             %                           when crit = Dopt and remove_inessential is true)

maxStepSearch = 10000; % Max number of loops for security
i = 0;

efficiency = 0;
while (i < maxStepSearch) && efficiency<eff
    i=i+1;
    Mold = Mnew;
    wold = wnew;
    % Step 1: find the direction
    maxDer = -Inf;
    maxDerIndex = 0;
    MoldInv = Mold\eye(p);  
    if crit=='Dopt'
        %----------------------------------------------------------------------
        % Remove inessential matrices ?
        if remove_inessential
            DER=zeros(1,K);
            for k = 1:K                     
                der = trAB(MoldInv, MM{k});
                DER(k)=der;
                if maxDer < der
                    maxDer = der;
                    maxDerIndex = k;
                end
            end
            epsilon=maxDer-p;
            threshold=p*(1+epsilon/2-sqrt(epsilon*(4+epsilon-4/p))/2 );
            i_keep=find(DER>=threshold);
            MM=MM(i_keep);
            indices=indices(i_keep);
            wold=wold(i_keep);
            K=length(i_keep);
            [~,ks] = max(DER(i_keep));
        else
            for k = 1:K                     
                der = trAB(MoldInv, MM{k});
                if maxDer < der
                    maxDer = der;
                    maxDerIndex = k;
                end
            end
            epsilon=maxDer-p;
            ks = maxDerIndex;    
        end
        efficiency=exp(-epsilon/p);
    elseif crit=='Lopt'
        % removing inessential matrices is more complicated and less
        % efficient than for D-optimal design and is not considered here
        for k = 1:K   
            KKtMm1=KK'*MoldInv;
            der = trace(KKtMm1*MM{k}*KKtMm1');
            if maxDer < der
                maxDer = der;
                maxDerIndex = k;
            end
        end
        ks = maxDerIndex;    
        efficiency=trace(KKtMm1*KK)/maxDer;
    end
    %----------------------------------------------------------------------
    
    % Step 2: search for best step by bisection
    Mks = MM{ks};
    wks = zeros(K, 1);
    wks(ks) = 1;
       
    a = 0;
    b = 1;
    
    for j = 1:maxBisec
        c = (a+b)/2;
        Mc=(1-c)*Mold + c*Mks;
       if crit=='Dopt' 
            Mcm1=inv(Mc);        
            further=trAB(Mcm1, (Mks-Mold)) >= 0;
        elseif crit=='Lopt' 
            KKtMcm1=KK'/Mc;
            further=trAB(KKtMcm1'*KKtMcm1, (Mks-Mold)) >= 0;
        end
        if further
            a = c;
        else 
            b = c;
        end
    end
    alpha = (a+b)/2; % Final step

    %alpha=1/(i+1);
      
    Mnew = (1-alpha) * Mold + alpha * Mks;
    wnew = (1-alpha) * wold + alpha * wks;
    
    if period<inf
        if rem(i,period)==0
            fprintf('\nEfficiency: %6.2f %%\n', 100*efficiency)
            fprintf([' \n'])
        end
    end    
end
Mopt=Mnew;
wopt=wnew;
    fprintf('\nEfficiency: %6.2f %%\n', 100*efficiency)
    fprintf('\nNb. of iterations: %6d\n', i)       
end