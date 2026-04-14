function [JML] = JML_stochastic(Y,controls,idx,nn,sig2,theta, Delta, Acell, B, C, D)
% function [JML] = JML_stochastic(Y,controls,idx,nn,sig2,theta, Delta, Acell, B, C, D)
% JML = - 2*log likelihood (up to a constant), for the stochastic model
% Ep_norm is just returned to check the quality (e.g., whiteness) of prediction errors

p=length(theta);
[n,~]=size(B);
[q,~]=size(C);
[qN,~]=size(Y);
N=qN/q;
A = 0*Acell{1};
for i = 1:p
    A = A+ theta(i) * Acell{i};
end
[F,G] = buildFG(A, B, Delta);
Sigma = lyap(A,D*D');
SigmaDelta = Sigma-F*Sigma*F';

Nexp=sum(nn);           % total number of experiments
Lnn=length(nn);

first_index_Y=1;
JML=0;

for i_input=1:Lnn
    nb_repeat=nn(i_input);  
    u = controls{idx(i_input)};
    for i_repeat=1:nb_repeat
        Yi=Y(:,first_index_Y); % observations for this input (q*N)*1
        %
        % Kalman filter 
        % --> Pk_given_km1 (N cells) and predk_given_km1 (N vectors q*1)
        [Pk_given_km1,predk_given_km1]=Kalman_1D(Yi,F,G,C,Sigma,zeros(n,1),u,SigmaDelta,sig2);
        for k=1:N
            Yk=Yi((k-1)*q+1:k*q);          % current q observations for t=t_k
            predk=predk_given_km1(:,k);    % current q predictions for t=t_k
            JML=JML+(Yk-predk)'*((Pk_given_km1{k}+sig2*eye(q))\(Yk-predk))+log(det(Pk_given_km1{k}+sig2*eye(q))); 
        end
        first_index_Y=first_index_Y+1;        
    end
end
JML=JML/qN/Nexp;
end