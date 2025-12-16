clear all;
addpath './Functions'
clc


rng(10)


%%%% Data

n = 2;
m = 1;
p = 2;
q = 1;

N=50;

Tf=2*pi;

prec = 500;   % must be even

Delta = Tf/N;


A0 = [-1, 1; 0, 1];
A1 = eye(n);
A2 = [0, -1; 1, 0];

Acol(:,:,1)=A1;
Acol(:,:,2)=A2;
Acol(:,:,3)=A0;

% sigma = 0.05*eye(n);
% Sigmap = 10^(-3);

sigma = 0.01*eye(n);
Sigmap = 10^(-4);


B = [1; 0];
C = [0, 1];

theta1 = -1/10;
theta2 = 2;

theta = [theta1,theta2];



v0 = ones(1,N);
v1 = cos([1:N]*2*pi/N);
v2 = cos(2*[1:N]*2*pi/N);
v3 = sin([1:N]*2*pi/N);
v4 = sin(2*[1:N]*2*pi/N);

density = 5; %% must be odd

controls = zeros(m,N,2);

counter=1;
for i0 = -density:2:density
    for i1 = -density:2:density
        for i2 = -density:2:density
            for i3 = -density:2:density
                for i4 = -density:2:density
                    u = (i0*v0+i1*v1+i2*v2+i3*v3+i4*v4);
                    u = u/sqrt(u*u'*N);
                    controls(:,:,counter) = u;
                    counter=counter+1;
                end
            end
        end
    end
end

[inputs,weights] = OptimalControls(Acol, B, C, theta, sigma, Sigmap, Delta, controls);


%%%% Creation of the p(p+1)/2 outputs



y = OutputGeneration(Acol, B, C, theta, sigma, Sigmap, Delta, inputs);

figure(1)
clf

X=[0:N]/N*Tf;

plot(X,[0,squeeze(y(:,:,1))],'b')
hold on
plot(X,[0,squeeze(y(:,:,2))],'r')
plot(X,[0,squeeze(y(:,:,3))],'g')

plot(Tf*[-0.05,1.05],[0,0],'k')
plot([0,0],1.05*[min(y,[],"all"),max(y,[],"all")],'k')



X0=X(1:end-1);
plot(X0,inputs(:,:,1),'b')
plot(X0,inputs(:,:,2),'r')
plot(X0,inputs(:,:,3),'g')


%%%% Differentiating the output


A=Acol(:,:,p+1);
for i = 1:p
    A=A+theta(i)*Acol(:,:,i);
end


F = buildF(A, Delta);
G = buildG(A, B, Delta, prec);
Sigma = buildSigma(A, sigma, Delta, prec);
SigmaBold = buildSigmaBold(C/2, C/2, F, Sigma, Sigmap, N);
SigmaBoldinv = inv(SigmaBold);


Bb = [B; zeros(n,m)];
Cb = [zeros(q,n), C];
Cbp = [C, zeros(q,n),];
sigmab = [sigma, zeros(n,n); zeros(n,n), zeros(n,n)];

display(G)

color=["--b";"--r";"--g"];


dJ = zeros(p,p*(p+1)/2);

Jaux = zeros(1,p*(p+1)/2);

for k = 1:p*(p+1)/2
    
    v = inputs(:,:,k);
    ybold = y(:,:,k);
    ybold = ybold(:);
    % eta(theta)
    
    eta = zeros(q,N);
    
    xt = zeros(n,1);
    for i=1:N
        xt = F*xt + G*v(:,i);
        eta(:,i) = C*xt;
    end
    etabold = eta(:);
    
    plot(X,[0;etabold],color(k))
  

    % Derivatives
    
    Jaux(k) = (ybold-etabold)'*SigmaBoldinv*(ybold-etabold)+log(det((SigmaBold)));
    

    for i = 1:p
        % Temporary elements
        Abi = [A, zeros(n,n); Acol(:,:,i), A];
        Fbi = buildF(Abi, Delta);
        Gbi = buildG(Abi, Bb, Delta, prec);
        Sigmabi = buildSigma(Abi, sigmab, Delta, prec);
    
        Sigmaboldbi = buildSigmaBold(Cb,Cbp,Fbi,Sigmabi,zeros(q,q),N);
    
        xbit = zeros(2*n,1);
        etab = zeros(q,N);
        for j= 1:N
            xbit = Fbi*xbit + Gbi*v(:,j);
            etab(:,j) = Cb*xbit;
        end
        etaboldbi = etab(:);
    
        dJ(i,k) = 2*(ybold-etabold)'*SigmaBoldinv*etaboldbi + (ybold-etabold)'*SigmaBoldinv*Sigmaboldbi*SigmaBoldinv*(ybold-etabold)+trace(SigmaBoldinv*Sigmaboldbi);
    
    end

end

J = Jaux * weights'
%gradJ = dJ * weights'


%[Jf, gradJf] = GLSDerivative(Acol, B, C, theta, sigma, Sigmap, Delta, inputs, weights, y)


Jv = GLSValue(Acol, B, C, theta, sigma, Sigmap, Delta, inputs, weights, y)

