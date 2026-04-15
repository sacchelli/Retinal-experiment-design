function [Pk_given_km1,predk_given_km1]=Kalman_1D(Y,F,G,C,P0,X0,u,Vnoise,sig2)
% function [Pk_given_km1,predk_given_km1]=Kalman_1D(Y,F,G,C,P0,X0,u,Vnoise,sig2)
% Kalman filter for the model X_{k+1}=F*X_k+G*u_k+W_k, W_k~N(0,Vnoise),
% and Y_k=C*X_k+e_k, e_k~N(0,sig2*eye(q))
% u (m*N) contains the N input vectors (dim = q) at times 1,...,N
% X_n is n*1 with X_0~N(X0,P0), G is n*m, C is q*n
% Pk_given_km1 = cell(N,1); each cell contains the posterior covariance of
% the prediction error
% predk_given_km1 is q*N, the kth column contains the prediction error at
% time k

%Epnorm=[];

[q,~]=size(C);
[~,N]=size(u);
predk_given_km1 = zeros(q,N);
Pk_given_km1 = cell(N,1);
% Initialisation k=0
Xk_k=X0;
Pk_k=P0;
for k=1:N
    Xk_km1=F*Xk_k+G*u(:,k);
    predk_given_km1(:,k)=C*Xk_km1;
    Pk_km1=F*Pk_k*F'+Vnoise;
    Pk_given_km1{k}=C*Pk_km1*C';
    % Kalman gain
    Kk=(Pk_km1*C')/(sig2*eye(q)+C*Pk_km1*C');
    Yk=Y((k-1)*q+1:k*q);
    Xk_k=Xk_km1+Kk*(Yk-predk_given_km1(:,k));
    Pk_k=Pk_km1-Kk*C*Pk_km1; 
    %Epnorm=[Epnorm (sqrtm(C*Pk_km1*C'+sig2*eye(q)))\(Yk-predk_given_km1(:,k))];
    %Pk_k=(Pk_k+Pk_k')/2;
end
% (eig(cov(Epnorm')))'
% pause
end