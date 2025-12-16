clear all;
addpath './Functions'
clc
warning('off', 'all')
%rng(10)
%rng(11)

%rng(315)
%rng(316)
rng(317)

%%%% Data

n = 3;
m = 1;
p = 3;
q = 1;

N = 20;

Tf=2*pi;

prec = 500;   % must be even

Delta = Tf/N;



% sigma = 0.5*eye(n);
% Sigmap = 0.5;

sigma = 0.1*eye(n);
Sigmap = 0.1;


% sigma = 0.01*eye(n);
% Sigmap = 10^(-4);


B = randn(n,1);
C = randn(1,n);

compteur = 0;
condition = false;
stop = 50;
while condition == false
    Acol = randn(n,n,p+1);
    thetaTrue = randn(1,p);
    A = Acol(:,:,p+1) + sum(Acol(:,:,1:p) .* reshape(thetaTrue, 1, 1, []), 3);
    condition = (max(real(eig(A)))<0) || (compteur > stop);
    compteur = compteur+1;
end
display(compteur)
display(eig(A))
%%


X=[1:N]/N*Tf;
v0 = ones(1,N);
v1 = cos(X/Tf*2*pi);
v2 = sin(X/Tf*2*pi);
v3 = cos(2*X/Tf*2*pi);
v4 = sin(2*X/Tf*2*pi);


density = 2; 

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
%%
%thetaInit = thetaTrue + 4 * abs(thetaTrue).*((rand(1,p))-0.5);

thetaInit = thetaTrue;

thetaOld = thetaInit;

%thetaOld = [ -415.0962 -533.3553]/1000;

thetaNew = thetaOld;

step=0.01;

estimate = [thetaInit];

estimateNaive = [thetaInit];

%%

loops = 1;

for k = 1:loops

    %if mod(k,5) == 1 
    [inputs,weights] = OptimalControls(Acol, B, C, thetaOld, sigma, Sigmap, Delta, controls);
    %end

    y = OutputGeneration(Acol, B, C, thetaTrue, sigma, Sigmap, Delta, inputs);
    
    %%%% Differentiating the output

    A=Acol(:,:,p+1);
    for i = 1:p
        A=A+thetaOld(i)*Acol(:,:,i);
    end


    F = buildF(A, Delta);
    G = buildG(A, B, Delta, prec);
    Sigma = buildSigma(A, sigma, Delta, prec);
    SigmaBold = buildSigmaBold(C/2, C, F, Sigma, Sigmap, N);


    [Jf, gradJf] = GLSDerivative(Acol, B, C, thetaOld, sigma, Sigmap, Delta, inputs, weights, y);

    
    %thetaNew = thetaOld + step * gradJf';
    normalizedGradJf = gradJf/(1000+norm(gradJf));
    %display(normalizedGradJf)
    thetaNew = thetaOld - step * normalizedGradJf';
    thetaOld = thetaNew;
    display(k/loops)
    estimate=[estimate;thetaOld];
end


% thetaOld = thetaInit;
% 
% %thetaOld = [ -415.0962 -533.3553]/1000;
% 
% thetaNew = thetaOld;
% 
% 
% 
% for k=1:loops
% 
%     inputs = controls(:,:,[fix(size(controls,3)/4),2*fix(size(controls,3)/4),3*fix(size(controls,3)/4)]);
% 
%     weights = [1/2,1/4,1/4];
% 
%     %[inputs,weights] = SubOptimalControl(Acol, B, C, thetaOld, sigma, Sigmap, Delta, controls);
% 
%     y = OutputGeneration(Acol, B, C, thetaTrue, sigma, Sigmap, Delta, inputs);
% 
% 
% 
%     %%%% Differentiating the output
% 
%     A=Acol(:,:,p+1);
%     for i = 1:p
%         A=A+thetaOld(i)*Acol(:,:,i);
%     end
% 
% 
%     F = buildF(A, Delta);
%     G = buildG(A, B, Delta, prec);
%     Sigma = buildSigma(A, sigma, Delta, prec);
%     SigmaBold = buildSigmaBold(C/2, C, F, Sigma, Sigmap, N);
% 
% 
%     [Jf, gradJf] = GLSDerivative(Acol, B, C, thetaOld, sigma, Sigmap, Delta, inputs, weights, y);
% 
% 
%     %thetaNew = thetaOld + step * gradJf';
%     normalizedGradJf = gradJf/(1000+norm(gradJf));
%     %display(normalizedGradJf)
%     thetaNew = thetaOld - step * normalizedGradJf';
%     thetaOld = thetaNew;
%     display(k/loops)
%     estimateNaive=[estimateNaive;thetaOld];
% end


% %%
% 
% % figure(1)
% % clf
% % plot(estimate(:,1),estimate(:,2),'b')
% % hold on
% % plot(estimateNaive(:,1),estimateNaive(:,2),'r')
% 
% 
% figure(1)
% clf
% 
% subplot(1,2,1)
% plot(sqrt(sum((estimate-ones(size(estimate,1),1)*thetaTrue).^2,2)),'b')
% hold on
% plot(sqrt(sum((estimateNaive-ones(size(estimate,1),1)*thetaTrue).^2,2)),'r')
% 
% subplot(1,2,2)
% semilogy(sqrt(sum((estimate-ones(size(estimate,1),1)*thetaTrue).^2,2)),'b')
% hold on
% semilogy(sqrt(sum((estimateNaive-ones(size(estimate,1),1)*thetaTrue).^2,2)),'r')
% %plot(sqrt((estimateNaive(:,1)-thetaTrue(1)).^2+(estimateNaive(:,2)-thetaTrue(2)).^2),'r')
% 


%%%%


T1min = thetaTrue(1)-abs(thetaTrue(1));
T1max = thetaTrue(1)+abs(thetaTrue(1));

T2min = thetaTrue(2)-abs(thetaTrue(2));
T2max = thetaTrue(2)+abs(thetaTrue(2));

nb=20;

dT1 = (T1max-T1min)/nb;
dT2 = (T2max-T2min)/nb;

T1 = T1min : dT1 : T1max;
T2 = T2min : dT2 : T2max;

[Theta1,Theta2] = meshgrid(T1,T2);
Jsurf=zeros(length(T2),length(T1));

[inputs, weights] = OptimalControlSingle(Acol, B, C, thetaTrue, sigma, Sigmap, Delta, controls);
y = OutputGeneration(Acol, B, C, thetaTrue, sigma, Sigmap, Delta, inputs);

tic
for i = 1:length(T1)
    for j = 1:length(T2)
        Jsurf(j,i) = GLSValue(Acol, B, C, [T1(i),T2(j), thetaTrue(3:end)], sigma, Sigmap, Delta, inputs, weights, y);
    end
end 
toc

inputs2 = controls(:,:,randi(size(controls,3),3,1));
y = OutputGeneration(Acol, B, C, thetaTrue, sigma, Sigmap, Delta, inputs2);

weights2 = rand(1,3);
weights2 = weights2/sum(weights2);

Jsurf2=zeros(length(T2),length(T1));
tic
for i = 1:length(T1)
    for j = 1:length(T2)
        Jsurf2(j,i) = GLSValue(Acol, B, C, [T1(i),T2(j), thetaTrue(3:end)], sigma, Sigmap, Delta, inputs2, weights2, y);
       
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
plot(thetaTrue(1),thetaTrue(2), 'r+', 'MarkerSize', 8, 'LineWidth', 2)

%hold on

subplot(1,2,2)
pcolor(Theta1, Theta2, log(Jsurf2-min(Jsurf2,[],"all")+1));
shading interp 
colorbar
hold on
plot(thetaTrue(1),thetaTrue(2), 'r+', 'MarkerSize', 8, 'LineWidth', 2)
%colormap( "parula");
%shading("interp");

