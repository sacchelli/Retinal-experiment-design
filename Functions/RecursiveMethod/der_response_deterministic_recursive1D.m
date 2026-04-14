function [dEtadti] = der_response_deterministic_recursive1D(itheta,controls,idx,nn,...
                                    theta, Delta, q, Acell, B, C, N)
% function [dEtadti] = der_response_deterministic_recursive1D(itheta,controls,idx,nn,...
%                                    theta, Delta, q, Acell, B, C, N)

% Computes the derivative of [Eta] = response_deterministic_recursive1D(controls,idx,nn,...
%                                   theta, Delta, prec, q, Acell, B, C, N) 
% with respect to the ith parameter theta(itheta)

if itheta<1 || itheta> length(theta)
    dEtadti=NaN;
    disp('itheta should be the index of a component of theta')
    return
end    

p=length(theta);
A = 0*Acell{1};
for i = 1:p
    A = A+ theta(i) * Acell{i};
end
[F,G] = buildFG(A, B, Delta);

[n,~]=size(B);
Ab = [A, zeros(n,n); Acell{itheta}, A];
Bb = [B; zeros(size(B))];
% Fb = buildF(Ab, Delta);
% Gb = buildG(Ab, Bb, Delta, prec);
[Fb,Gb] = buildFG(Ab,Bb, Delta);
dFdti=Fb(n+1:2*n,1:n); % derivative of F with respect to theta(itheta)
    % % check derivative
    % dFdti_approx=1e9*(buildF(A+1e-9*Acell{itheta}, Delta) - buildF(A, Delta));
    % error_F=max(max(abs(dFdti-dFdti_approx)))
    % % pause
dGdti=Gb(n+1:2*n,:); % derivative of G with respect to theta(itheta)
    % % check derivative
    % dGdti_approx=1e9*(buildG(A+1e-9*Acell{itheta}, B, Delta, prec) - buildG(A, B, Delta, prec));
    % error_G=max(max(abs(dGdti-dGdti_approx)))
    % pause
Lnn=length(nn);
dEtadti=zeros(q*N,Lnn);    % derivatives of vectors of model responses
for i_input=1:Lnn
    % compute eta (deterministic response, length q*N) for i_input
    % together with its derivative with respect to theta(itheta)
    u = controls{idx(i_input)};
    Xk=G*u(:,1);
    %dXkdti=zeros(size(Xk));
    dXkdti=dGdti*u(:,1);
    dEtadti(1:q,i_input)=C*dXkdti;
    for k=2:N
        dXkdti=F*dXkdti+dFdti*Xk+dGdti*u(:,k);
        Xk=F*Xk+G*u(:,k);
        dEtadti((k-1)*q+1:k*q,i_input)=C*dXkdti;
    end
end 

end