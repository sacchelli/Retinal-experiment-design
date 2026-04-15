function [Y] = data_stochastic_1D_rec(controls,idx,nn,sig2,...
                                    theta, Delta, Acell, B, C, D, N)
% function [Y] = data_stochastic_1D_rec(controls,idx,nn,sig2,...
%                                    theta, Delta, Acell, B, C, D, N)
% for the model X_{k+1}=F*X_k+G*u_k+W_k, W_k~N(0,SigmaDelta), 
% with X_0 in stationary condition, X_0~N(0,Sigma) 
% and Y_k=C*X_k+e_k, e_k~N(0,sig2*eye(q))
% F, G, Sigma and SigmaDelta computed by discretization at time-steps Delat 
% of the time continous model % dX(t) = A*X(t)+B*u(t)+D dW(t), 
% with W a standard Brownian motion. 
% idx = indices of inputs (elements of controls) to be used, with
% nn_i>0 = nb. of repetitions of input = controls{idx(i)}
% Returns a matrix of observations Y = matrix (q*N)*Nexp, with 
% Y(:,i) the vector of observations at the ith input 
% (taken in order of appearance in nn, with repetitions), 
% and Nexp = sum(nn) the number of experiments for parameters theta

p=length(theta);
[n,~]=size(B);
[q,~]=size(C);
A = 0*Acell{1};
for i = 1:p
    A = A+ theta(i) * Acell{i};
end
[F,G] = buildFG(A, B, Delta);


Sigma = lyap(A,D*D');
SigmaDelta = Sigma-F*Sigma*F';

qN=q*N;
Nexp=sum(nn);           % total number of experiments
Lnn=length(nn);
qN
Nexp
Y=zeros(qN,Nexp);  % matrix of Nexp vectors of (q*N)*1 observations
first_index_Y=1;
for i_input=1:Lnn
    nb_repeat=nn(i_input);   
    % compute eta (deterministic response, length q*N) for i_input
    u = controls{idx(i_input)};
    for i_repeat=1:nb_repeat
        % generate observations 
        Xk=mvnrnd(zeros(n,1), Sigma, 1); % random initial state
        Xk=Xk';
        for k=1:N
            Wk=mvnrnd(zeros(n,1), SigmaDelta, 1);
            Xk=F*Xk+G*u(:,k)+Wk';
            Y((k-1)*q+1:k*q,first_index_Y)=C*Xk+sqrt(sig2)*randn(q,1);
        end
        first_index_Y=first_index_Y+1;
    end
end
end