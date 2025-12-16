clc

n=50;
q=1;

N=20;


F = rand(n,n)/n;

C = rand(q,n);

Sigma = rand(n,n);
Sigma = Sigma*Sigma';
Sigmap = .1*eye(q);




tic
M0 = buildSigmaBold(C, C/2, F, Sigma, Sigmap, N);
toc

tic
M1 = buildSigmaBoldCutOff(C, C/2, F, Sigma, Sigmap, N);
toc



%max(abs(M1-M2)./abs(M1),[],"all")

%display(max(abs(M0-M1)./abs(M0),[],"all"))

%display(max(abs(M0-M2)./abs(M0),[],"all"))


%abs(M0-M1)






%%

sdn=10;

dx = 1/sdn;
dy = 1/sdn;

Tf=2*pi;

N=100;
Delta = Tf/N;

thetaU = pi/3;
lambdaU = 1;
vU = 5;

figure(1)
clf


input=zeros(sdn^2,N);

compteur = 0;

for k = 1:N
    compteur = compteur+1;
    temp=zeros(sdn,sdn);
    for i = 1:sdn
        for j = 1:sdn
            temp(i,j)=1+cos(2*pi*i*dx*cos(thetaU)+2*pi*j*dy*sin(thetaU)+k*Delta*vU);
        end
    end
    input(:,k) = temp(:);
    imagesc(temp)
    drawnow
    pause(0.05)
end