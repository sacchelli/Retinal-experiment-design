function [Eta] = response_deterministic_recursive1D(controls,idx,nn,...
                                    theta, Delta, q, Acell, B, C, N)
% function [Eta] = response_deterministic_recursive1D(controls,idx,nn,...
%                                     theta, Delta, prec, q, Acell, B, C, N)

p=length(theta);
A = 0*Acell{1};
for i = 1:p
    A = A+ theta(i) * Acell{i};
end
[F,G] = buildFG(A, B, Delta);

Lnn=length(nn);
Eta=zeros(q*N,Lnn);    % vectors of model responses
for i_input=1:Lnn
    % compute eta (deterministic response, length q*N) for i_input
    u = controls{idx(i_input)};
    Xk=G*u(:,1);
    Eta(1:q,i_input)=C*Xk;
    for k=2:N
        Xk=F*Xk+G*u(:,k);
        Eta((k-1)*q+1:k*q,i_input)=C*Xk;
    end
end 
end