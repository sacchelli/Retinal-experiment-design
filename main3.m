clear all;
addpath './Functions'
clc


%rng(10)
%rng(11)

%rng(201)
rng(202)

%%%% Data

n = 2;
m = 1;
p = 2;
q = 1;

N=20;

Tf=2*pi;

prec = 500;   % must be even

Delta = Tf/N;


A0 = [-1, 1; 0, 1];
A1 = eye(n);
A2 = [0, -1; 1, 0];

Acol(:,:,1)=A1;
Acol(:,:,2)=A2;
Acol(:,:,3)=A0;

% sigma = 0.5*eye(n);
% Sigmap = 0.05;

sigma = 0.01*eye(n);
Sigmap = 0.1;

% sigma = 0.01*eye(n);
% Sigmap = 10^(-4);


B = [1; 0];
C = [0, 1];

% theta1 = -1/10;
theta1 = -1;
theta2 = 2;

theta = [theta1,theta2];


% X=[1:N]/N*Tf;
% v0 = ones(1,N);
% v1 = cos(X/Tf*2*pi);
% v2 = sin(X/Tf*2*pi);
% v3 = cos(2*X/Tf*2*pi);
% v4 = sin(2*X/Tf*2*pi);
% v5 = cos(3*X/Tf*2*pi);
% v6 = sin(3*X/Tf*2*pi);
% 
% 
% density = 2; %% must be odd
% 
% controls = zeros(m,N,2);
% 
% counter = 1;
% for i0 = -density:1:density
%     for i1 = -density:1:density
%         for i2 = -density:1:density
%             for i3 = -density:1:density
%                 for i4 = -density:1:density
%                     for i5 = -density:1:density
%                         for i6 = -density:1:density
%                             u = (i0*v0+i1*v1+i2*v2+i3*v3+i4*v4+i5*v5+i6*v6);
% 
%                             L2u = max(u*u'/N,0.000001);
%                             u = u/sqrt(L2u);
%                             controls(:,:,counter) = u;
%                             counter = counter+1;
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end


X=[1:N]/N*Tf;
v0 = ones(1,N);
v1 = cos(X/Tf*2*pi);
v2 = sin(X/Tf*2*pi);
v3 = cos(2*X/Tf*2*pi);
v4 = sin(2*X/Tf*2*pi);


density = 3; 

controls = zeros(m,N,2);

counter = 1;
for i0 = -density:1:density
    for i1 = -density:1:density
        for i2 = -density:1:density
            for i3 = -density:1:density
                for i4 = -density:1:density
                            u = (i0*v0+i1*v1+i2*v2+i3*v3+i4*v4);

                            L2u = max(u*u'/N,0.000001);
                            u = u/sqrt(L2u);
                            controls(:,:,counter) = u;
                            counter = counter+1;
                end
            end
        end
    end
end

[inputs,weights] = OptimalControls(Acol, B, C, theta, sigma, Sigmap, Delta, controls);


%inputs = controls(:,:,[randi(size(controls,3)),randi(size(controls,3)),randi(size(controls,3))])

%%%% Creation of the p(p+1)/2 outputs


y = OutputGeneration(Acol, B, C, theta, sigma, Sigmap, Delta, inputs);

figure(1)
clf

colors = [0.20 0.65 0.80;   % blue
          0.90 0.60 0.00;   % orange
          0.7 0.3 0.90];  % purple
colors = ["#1E88E5";"#FFC107";"#D81B60"]

X=[0:N]/N*Tf;

plot(X,[0,squeeze(y(:,:,1))],"Color",colors(1,:),"LineWidth",1.)
hold on
plot(X,[0,squeeze(y(:,:,2))],"Color",colors(2,:),"LineWidth",1.)
plot(X,[0,squeeze(y(:,:,3))],"Color",colors(3,:),"LineWidth",1.2)


xlim([-0.02*Tf,1.02*Tf])
xline(0)
yline(0)



X0=X(1:end-1);
plot(X0,inputs(:,:,1),'--',"Color",colors(1,:),"LineWidth",1.)
plot(X0,inputs(:,:,2),'--',"Color",colors(2,:),"LineWidth",1.)
plot(X0,inputs(:,:,3),'--',"Color",colors(3,:),"LineWidth",1.2)


%%%% Differentiating the output

A=Acol(:,:,p+1);
for i = 1:p
    A=A+theta(i)*Acol(:,:,i);
end


F = buildF(A, Delta);
G = buildG(A, B, Delta, prec);
Sigma = buildSigma(A, sigma, Delta, prec);
SigmaBold = buildSigmaBold(C/2, C, F, Sigma, Sigmap, N);
%SigmaBoldinv = inv(SigmaBold);





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
    
    plot(X,[0;etabold],"Color",colors(k,:),"LineWidth",1.)
end

%profile on
tic

[Jf, gradJf] = GLSDerivative(Acol, B, C, theta, sigma, Sigmap, Delta, inputs, weights, y)

toc

%profile viewer

tic
Jv = GLSValue(Acol, B, C, theta, sigma, Sigmap, Delta, inputs, weights, y)
toc


%% 

%%%% Plotting


% T1min = theta1-30*abs(theta1);
% T1max = theta1+30*abs(theta1);
% 
% T2min = theta2-0.4*abs(theta2);
% T2max = theta2+0.4*abs(theta2);

T1min = -2;
T1max = -0.0001;

T2min = 0;
T2max = 3;

nb=60;

dT1 = (T1max-T1min)/nb;
dT2 = (T2max-T2min)/nb;

T1 = T1min : dT1 : T1max;
T2 = T2min : dT2 : T2max;

[Theta1,Theta2] = meshgrid(T1,T2);

Jsurf=zeros(length(T2),length(T1));
tic
for i = 1:length(T1)
    for j = 1:length(T2)
        Abis = A0+ T1(i)*A1+T2(j)*A2;       
        %if max(real(eig(Abis)))>0
        %    Jsurf(j,i)=Inf;
        %else
            Jsurf(j,i) = GLSValue(Acol, B, C, [T1(i),T2(j)], sigma, Sigmap, Delta, inputs, weights, y);
        %end
    end
end 
toc





inputs2 = controls(:,:,randi(size(controls,3),3,1));

weights2 = rand(1,3);
weights2 = weights2/sum(weights2);

Jsurf2=zeros(length(T2),length(T1));
tic
for i = 1:length(T1)
    for j = 1:length(T2)
        Abis = A0+ T1(i)*A1+T2(j)*A2;       
        %if max(real(eig(Abis)))>0
        %    Jsurf(j,i)=Inf;
        %else
            Jsurf2(j,i) = GLSValue(Acol, B, C, [T1(i),T2(j)], sigma, Sigmap, Delta, inputs2, weights2, y);
        %end
    end
end 
toc



figure(2)
clf
subplot(1,2,1)

surf(Theta1, Theta2, Jsurf);
colormap("parula");
shading("interp");

%hold on

subplot(1,2,2)
surf(Theta1, Theta2, Jsurf2);
colormap( "parula");
shading("interp");

%%


figure(3)
clf
subplot(1,2,1)

pcolor(Theta1, Theta2, log(Jsurf-min(Jsurf,[],"all")+1));
shading interp 
colorbar
hold on
plot(theta1, theta2, 'r+', 'MarkerSize', 8, 'LineWidth', 2)

%hold on

subplot(1,2,2)
pcolor(Theta1, Theta2, log(Jsurf2-min(Jsurf2,[],"all")+1));
shading interp 
colorbar
hold on
plot(theta1, theta2, 'r+', 'MarkerSize', 8, 'LineWidth', 2)
%colormap( "parula");
%shading("interp");

