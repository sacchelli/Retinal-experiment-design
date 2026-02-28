function [Eta] = response_deterministic_1D(controls,idx,nn,...
                                    theta, Delta, prec, q, Acell, B, C, N)
% function [Eta] = response_deterministic_1Dcontrols,idx,nn,...
%                                     theta, Delta, prec, q, Acell, B, C, N)

p=length(theta);
A = 0*Acell{1};
for i = 1:p
    A = A+ theta(i) * Acell{i};
end
[F,G] = buildFG(A,B, Delta);
H = buildH(C, F, G, N);

Lnn=length(nn);
Eta=zeros(q*N,Lnn);    % vectors of model responses
for i_input=1:Lnn
    % compute eta (deterministic response, length q*N) for i_input
    u = controls{idx(i_input)};
    u = u(:);
    Eta(:,i_input)=H*u;
end 
end