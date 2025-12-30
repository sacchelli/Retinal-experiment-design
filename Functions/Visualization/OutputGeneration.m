function outputs = OutputGeneration(Acol, B, C, theta, sigma, Sigmap, Delta, controls)

% I assume the following shapes: 
% Acol of shape n x n x (p+1) with A_0 = Acol(:,:,p+1), A_i = Acol(:,:,i)
% theta of shape 1 x p
% B of shape n x m
% C of shape q x n
% sigma, Sigma of shape n x n
% Sigmap of shape q x q
% controls of shape m x N x nU
%
% outputs should be of shape q x N x nU

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
CholSigma = chol(Sigma, 'lower');    % Cholesky factor, Sigma = L*L'
CholSigmap = chol(Sigmap, 'lower');

%%%% generation of the stochatic trajectories
ybold = zeros(q, N, nU);

for k = 1:nU
    v = controls(:,:,k);
    xt = zeros(n,1);
    yAux = zeros(q,N);
    for i=1:N
        Wk = CholSigma*randn(n, 1);
        epsk = CholSigmap*randn(q, 1);
        xt = F*xt + G*v(:, i) + Wk;
        yAux(:,i) = C*xt + epsk;
    end    
    ybold(:,:,k) = yAux; 
end



outputs = ybold;

end