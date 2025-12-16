function [J, gradJ] = GLSValue(Acol, B, C, theta, sigma, Sigmap, Delta, controls, weights, y)

% I assume the following shapes: 
% Acol of shape n x n x (p+1) with A_0 = Acol(:,:,p+1), A_i = Acol(:,:,i)
% theta of shape 1 x p
% B of shape n x m
% C of shape q x n
% sigma, Sigma of shape n x n
% Sigmap of shape q x q
% controls of shape m x N x nU
% weights of shape 1 x nU
% y of shape q x N x nU

n = size(B,1);
p = length(theta);
q = size(C,1);
N = size(controls, 2);
nU = size(controls, 3);

prec = 500;   % parameter for the integration by quadrature, must be even


A=Acol(:,:,p+1);
for i = 1:p
    A=A+theta(i)*Acol(:,:,i);
end


%%%% Construction of elements

F = buildF(A, Delta);
G = buildG(A, B, Delta, prec);
Sigma = buildSigma(A, sigma, Delta, prec);
SigmaBold = buildSigmaBold(C/2, C/2, F, Sigma, Sigmap, N);
%SigmaBoldinv = inv(SigmaBold);

Jaux = zeros(1,nU);


for k = 1:nU
    
    v = controls(:,:,k);

    ybold = y(:,:,k);
    ybold=ybold(:);
    
    

    % eta(theta)
    
    eta = zeros(q,N);
    
    xt = zeros(n,1);
    for i=1:N
        xt = F*xt + G*v(:,i);
        eta(:,i) = C*xt;
    end
    
    etabold = eta(:);

    % Derivatives
    
    Jaux(k) = real((ybold-etabold)'*(SigmaBold\(ybold-etabold))+logdet(SigmaBold));
    

end

J = Jaux * weights';



end