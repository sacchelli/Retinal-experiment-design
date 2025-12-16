clear all;
addpath './Functions'
clc

tic

%%%% Data

n = 2;
m = 1;
p = 2;
q = 1;

N = 10;

prec = 500;   % must be even

Delta = 1;


A0 = [-1, 1; 0, 1];
A1 = eye(n);
A2 = [0, -1; 1, 0];

sigma = 0.01*eye(n);
Sigmap = 10^(-4);


B = [1; 0];
C = [0, 1];

% theta1 = -1/10;
theta1 = -1;
theta2 = 2;


theta1 = -0.4151;
theta2 = -0.5334;

A = A0+theta1*A1+theta2*A2;


%%%% Construction of elements


F = buildF(A, Delta);
G = buildG(A, B, Delta,prec);
Sigma = buildSigma(A, sigma, Delta, prec);

Ab1 = [A, zeros(n,n); A1, A];
Ab2 = [A, zeros(n,n); A2, A];

Bb = [B; zeros(n,m)];
Cb = [zeros(q,n), C];
Cbp = [C, zeros(q,n),];
sigmab = [sigma, zeros(n,n); zeros(n,n), zeros(n,n)];

Fb1 = buildF(Ab1, Delta);
Fb2 = buildF(Ab2, Delta);

Gb1 = buildG(Ab1, Bb, Delta, prec);
Gb2 = buildG(Ab2, Bb, Delta, prec);

Sigmab1 = buildSigma(Ab1, sigmab, Delta, prec);
Sigmab2 = buildSigma(Ab2, sigmab, Delta, prec);

Hb1 = buildH(Cb, Fb1, Gb1, N);
Hb2 = buildH(Cb, Fb2, Gb2, N);

SigmaBold = buildSigmaBold(C,C/2,F,Sigma,Sigmap,N);

SigmaBoldinv = inv(SigmaBold);

dS1 = buildSigmaBold(Cb,Cbp,Fb1,Sigmab1,zeros(q,q),N);
dS2 = buildSigmaBold(Cb,Cbp,Fb2,Sigmab1,zeros(q,q),N);

Q11 = 0.5*trace(SigmaBoldinv*dS1*SigmaBoldinv*dS1);
Q22 = 0.5*trace(SigmaBoldinv*dS2*SigmaBoldinv*dS2);
Q12 = 0.5*trace(SigmaBoldinv*dS1*SigmaBoldinv*dS2);

L11 = Hb1'*SigmaBoldinv*Hb1;
L22 = Hb2'*SigmaBoldinv*Hb2;
L12 = (Hb1'*SigmaBoldinv*Hb2+Hb2'*SigmaBoldinv*Hb1)/2;





%%%% Collection of controls

%%% Test family 1

% controls = [0;1];
% for i = 1:N-1
%     controls = [zeros(size(controls,1),1), controls; ones(size(controls,1),1), controls];
% end
% 
% controls=controls(2:end,:);
% 
% for i = size(controls,1)
%     controls(i,:) = controls(i,:)/sqrt(controls(i,:)*controls(i,:)');
% end


%%% Test family 2

v0 = ones(1,N);
v1 = cos([1:N]*2*pi/N);
v2 = cos(2*[1:N]*2*pi/N);
v3 = sin([1:N]*2*pi/N);
v4 = sin(2*[1:N]*2*pi/N);

density = 3; %% must be odd

controls = [];
for i0 = -density:2:density
    for i1 = -density:2:density
        for i2 = -density:2:density
            for i3 = -density:2:density
                for i4 = -density:2:density
                    u = (i0*v0+i1*v1+i2*v2+i3*v3+i4*v4);
                    L2u = max(u*u'/N,0.000001);
                    %u = u/sqrt(u*u');
                    u = u/sqrt(L2u);
                    controls=[controls;u];
                end
            end
        end
    end
end

%%% Test family 3

% controls = [];
% 
% for i0 = 1:100
%     u = cos([1:N]/N*log(1+i0));
%     u = u/sqrt(u*u');
%     controls=[controls;u];
% end


%%%% Collection of Fisher matrices

M = zeros(p,p,size(controls,1));

for i = 1:size(controls,1)
    u=controls(i,:)';
    M(1,1,i) = u'*L11*u+Q11;
    M(1,2,i) = u'*L12*u+Q12;
    M(2,2,i) = u'*L22*u+Q22;
    M(2,1,i) = M(1,2,i);
end

% Idea for upper matrix extraction:

% mask = triu(true(p));                 % logical mask for upper triangular
% rowVec = reshape(M(i,:,:), p, p);     % turn slice into p×p
% rowVec = rowVec(mask).';   

FisherCollection=[squeeze(M(1,1,:)),squeeze(M(1,2,:)),squeeze(M(2,2,:))];




ConvexHullGraph = convhulln(FisherCollection);

vertices = unique(ConvexHullGraph(:))';

%%%% Using a priori estimates on the score

maxscorevertex=-Inf;
maxscoreverpos=0;

for k = vertices
    if maxscorevertex<log(det(M(:,:,k)))
        maxscorevertex = log(det(M(:,:,k)));
        maxscoreverpos=k;
    end
    %maxscorevertex = max(maxscorevertex,log(det(M(:,:,k))));

end


display(maxscorevertex)
display(maxscoreverpos)







apriorifaces = [];

for i = 1:size(ConvexHullGraph,1)

    score= UpperBoundFace(M(:,:,ConvexHullGraph(i,:)));
    
    if score >=  maxscorevertex
        apriorifaces=[apriorifaces,i];
    end
end

ConvHullCap=ConvexHullGraph(apriorifaces,:)



%%%% Computing the score of faces on the cap

scoreMax = -Inf;

faceMax = 0;

wMax = 0;

for i = 1:size(ConvHullCap,1)

    [score, weights] = FaceScoreDual(M(:,:,ConvHullCap(i,:)));
    
    if score > scoreMax
        scoreMax = score;
        faceMax = ConvHullCap(i,:);
        wMax = weights;
    end
    

end

display(scoreMax)
display(faceMax)



%%


%%%% Display of results

ConeCoordinates = [squeeze(M(1,2,:)), squeeze(M(1,1,:)-M(2,2,:))/2, squeeze(M(1,1,:)+M(2,2,:))/2];

coordMax=max(max(abs(ConeCoordinates(:,1:2))));

figure(1)
clf()

[x, y] = meshgrid(-coordMax:coordMax/50:coordMax, -coordMax:coordMax/50:coordMax);

% Compute z for upper cone
z_pos = sqrt(x.^2 + y.^2); 

% Plot both surfaces

surf(x, y, z_pos, 'FaceAlpha', 0.7); 
colormap("parula")
shading interp
hold on


trisurf(ConvexHullGraph,ConeCoordinates(:,1),ConeCoordinates(:,2),ConeCoordinates(:,3),'FaceColor','yellow')
hold on

trisurf(ConvHullCap,ConeCoordinates(:,1),ConeCoordinates(:,2),ConeCoordinates(:,3),'FaceColor','cyan')

%scatter3(M(:,1,1),M(:,1,2),M(:,2,2))

figure(2)
clf



coordMaxCapX=max(ConeCoordinates(unique(ConvHullCap(:)),1));
coordMinCapX=min(ConeCoordinates(unique(ConvHullCap(:)),1));


coordMaxCapY=max(ConeCoordinates(unique(ConvHullCap(:)),2));
coordMinCapY=min(ConeCoordinates(unique(ConvHullCap(:)),2));

[xbis, ybis] = meshgrid(coordMinCapX:(coordMaxCapX-coordMinCapX)/50:coordMaxCapX, coordMinCapY:(coordMaxCapY-coordMinCapY)/50:coordMaxCapY);

z_posbis = sqrt(xbis.^2 + ybis.^2+exp(maxscorevertex));

surf(xbis, ybis, z_posbis, 'FaceAlpha', 0.7); 
colormap("parula")
shading interp
hold on 

trisurf(ConvHullCap,ConeCoordinates(:,1),ConeCoordinates(:,2),ConeCoordinates(:,3),'FaceColor','cyan')


toc