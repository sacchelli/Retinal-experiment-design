function [Y] = data_deterministic_1D_rec(controls,idx,nn,sig2,...
                                    theta, Delta, Acell, B, C, N)
% function [Y] = data_deterministic_1D_rec(controls,idx,nn,sig2,...
%                                     theta, Delta, Acell, B, C, N)
% for the model X_{k+1} = F*X_k+G*u_k, with X_0 = 0
% and Y_k = C*X_k+e_k, e_k~N(0,sig2*eye(q))
% idx = indices of inputs (elements of controls) to be used, with
% nn_i>0 = nb. of repetitions of input = controls{idx(i)}
% Returns a matrix of observations Y = matrix (q*N)*Nexp, with 
% Returns a matrix of observations Y = matrix (q*N)*Nexp, with 
% Y(:,i) the vector of observations at the ith input 
% (taken in order of appearance in nn, with repetitions), 
% and Nexp = sum(nn) the number of experiments for parameters theta

p=length(theta);
[q,~]=size(C);
A = 0*Acell{1};
for i = 1:p
    A = A+ theta(i) * Acell{i};
end
[F,G] = buildFG(A, B, Delta);

Lnn=length(nn);
Nexp=sum(nn);
Y=zeros(q*N,Nexp);  % matrix of Nexp vectors of (q*N)*1 observations
first_index_Y=1;
for i_input=1:Lnn
    nb_repeat=nn(i_input);   
    % compute eta (deterministic response, length q*N) for i_input
    u = controls{idx(i_input)};
    Xk=G*u(:,1);
    Eta=zeros(q*N,1);
    Eta(1:q)=C*Xk;
    for k=2:N
        Xk=F*Xk+G*u(:,k);
        Eta((k-1)*q+1:k*q)=C*Xk;
    end
    for i_repeat=1:nb_repeat
        % generate observations = eta + noise normal(0,sig2)
        Y(:,first_index_Y)=Eta+sqrt(sig2)*randn(N*q,1); 
        first_index_Y=first_index_Y+1;
    end
end 
end