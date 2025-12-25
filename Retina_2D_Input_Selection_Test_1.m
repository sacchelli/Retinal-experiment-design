addpath './Functions'
clc

timeTakenStart = tic;

% Pick the rng seed if necessary.

seed = 1;

rng(seed);

% Should we use cutoff in the large, near diagonal, matrices ?

cutoff = true;

% if yes, then which depth ? It will be depth * q. 20 is a good choice
% under current data.

depth = 20;



%%%% Data // following the reference.


retinalWidth = 2; % mm 

%delta = 0.00389; % mm
sigmaBG = 0.15; % mm
sigmaAG = 0.15; % mm

tauB = 57.5139; % ms
tauA = 86.6298; % ms
tauG = 54.8947; % ms

wm = 0.0128109; % kHz
wp = 0.0188413; % kHz
wBG = 0.0131284; % kHz 
wAG = 0.0137853; % kHz 

%%%% Numerical data

sdn = 10; % size of the edge of one layer

dn = sdn^2; % size of one layer : needs to be a square

n = 3 * dn;

m = dn; % number of inputs

p = 7; % number of parameters

q = dn; % number of outputs

N = 100; % Number of measurements


delta = retinalWidth/dn; % distance between cells


%%%% Model


Gamma = adjacencyLaplacian2DGraph(sdn);
% Gamma = absoluteLaplacian2DGraph(sdn);
% Gamma = laplacian2DGraph(sdn);
% Gamma = closestNeighbor2DGraph(sdn);

GammaBG = gaussianPooling2DGraph(sdn,delta,sigmaBG);
GammaAG = gaussianPooling2DGraph(sdn,delta,sigmaAG);


id = speye(dn);
zz = sparse(dn,dn);


Acol = cell(p+1, 1);

Acol{1} = [-id , zz, zz ; zz, zz, zz ; zz, zz, zz];
Acol{2} = [zz , zz, zz ; zz, -id, zz ; zz, zz, zz];
Acol{3} = [zz , zz, zz ; zz, zz, zz ; zz, zz, -id];

Acol{4} = [zz , -Gamma, zz ; zz, zz, zz ; zz, zz, zz];
Acol{5} = [zz , zz, zz ; Gamma, zz, zz ; zz, zz, zz];

Acol{6} = [zz , zz, zz ; zz, zz, zz ; GammaBG, zz, zz];
Acol{7} = [zz , zz, zz ; zz, zz, zz ; zz, GammaAG, zz];

Acol{8} = 1/tauB*Acol{1} + 1/tauA*Acol{2} + 1/tauG*Acol{3} ...
    + wm*Acol{4} + wp*Acol{5} + wBG*Acol{6} + wAG*Acol{7};



Delta = max([tauA,tauB,tauG])/5;

Tf = N * Delta; % Time of the experiment

prec = 500;   % must be even



%%%% Determination of 

sigma = 0.1*eye(n);
Sigmap = 0.1*eye(q);

B = sparse([eye(dn,m);zeros(2*dn,m)]);

C = sparse([zeros(q,2*dn),eye(q,dn)]);


%%% Contraction of A

thetaData = [1/tauB, 1/tauA, 1/tauG, wm, wp, wBG, wAG];

%thetaTrue = (rand(1,p)-1/2)/2.*thetaData;
%thetaTrue = 0.*thetaData;
thetaTrue = thetaData;

A = Acol{p+1};
for i = 1:p
    A = A + thetaTrue(i) * Acol{i};
end

maxReEigA = max(real(eigs(A)));



fprintf(['Data created, Max Re(λ(A)) = ', num2str(maxReEigA), '. \n'])


%%%% Construction of elements

fprintf(['\n','Dynamics and variance pre-computation... '])

tStepStart = tic;

F = buildF(A, Delta);
G = buildGASparse(A, B, Delta, prec);
Sigma = buildSigmaASparse(A, sigma, Delta, prec);
SigmaBold = buildSigmaBold(C, C/2, F, Sigma, Sigmap, N);
SigmaBoldinv = inv(SigmaBold);


% Use this to accelerate the computations by killing off diagonal terms:

if cutoff
    SigmaBoldinv = sparse(bandDiagonalFilter(SigmaBoldinv,depth*q)); 
end


tStepEnd = toc(tStepStart);

fprintf(['done (', num2str(tStepEnd,3) ,'s). \n'])

fprintf(['\n','Parameter sensitivities pre-computation... \n'])

%
Hb = cell(p, 1);
dSb = cell(p, 1);

for i = 1:p
    fprintf(['Parameter ',num2str(i),'/',num2str(p),' '])
    tStepStart = tic;
    
    Ab = [A, zeros(n,n); Acol{i}, A];
    
    Bb = [B; zeros(n,m)];
    Cb = [zeros(q,n), C];
    Cbp = [C, zeros(q,n)];
    sigmab = [sigma, zeros(n,n); zeros(n,n), zeros(n,n)];

    Fb = buildF(Ab, Delta);
    
    Gb = buildGASparse(Ab, Bb, Delta, prec);
    
    Sigmab = buildSigmaASparse(Ab, sigmab, Delta, prec);
    
    Hb{i} = buildH(Cb, Fb, Gb, N);
    dSb{i} = buildSigmaBold(Cb, Cbp, Fb, Sigmab, zeros(q,q), N);
    if cutoff
        Hb{i} = sparse(bandDiagonalFilter(Hb{i},depth*q));
        dSb{i} = sparse(bandDiagonalFilter(dSb{i},depth*q));
    end


    tStepEnd = toc(tStepStart);

    fprintf(['done (', num2str(tStepEnd,3) ,'s). \n'])
end

timeTakenCheck = toc(timeTakenStart);

fprintf(['\n','Total computation time up to now: ', num2str(timeTakenCheck,3) ,'s\n'])

fprintf(['\n','Q matrix computation... \n'])


Q = zeros(p,p);
halfQ = cell(p, 1);

for i = 1:p
    fprintf(['Parameter ',num2str(i),'/',num2str(p)])
    tStepStart = tic;
    halfQ{i} = full(SigmaBoldinv) * full(dSb{i});
    tStepEnd = toc(tStepStart);
    fprintf([' done (', num2str(tStepEnd,3) ,'s).\n'])
end

for i = 1:p
    for j = 1:i
        Q(i,j) = 0.5*trAB(halfQ{i}, halfQ{j});
        Q(j,i) = Q(i,j);
    end
end


timeTakenCheck = toc(timeTakenStart);

fprintf(['\n','Total computation time up to now: ', num2str(timeTakenCheck,3) ,'s\n'])

%%
%%%% Collection of controls

%%% Test family 

fprintf(['\n','Input collection... '])

tStepStart = tic;


cMax = 2/Tf;   % max wave velocity
cMin = 0; % min wave velocity
kMax = 2*pi/retinalWidth*(dn-1)/2;   % max wave number
kMin = 2*pi/(retinalWidth);   % min wave number

cN = 10;
kN = 10;
thetaN = 10;


[controls2DWave, inputParameters2Dwave] = generate2DPlaneWaveControls(sdn, N, Delta, retinalWidth, thetaN, cMin, cMax, cN, kMin, kMax, kN);

[controlsRandom, inputParametersRandom] = generate2DRandomControls(sdn, N, 100, seed);

% constant input

%counter = length(controls2DWave) + 1;

counter = length(controls2DWave) + length(controlsRandom) + 1;

controls = cell(counter,1);

%inputParameters = [inputParameters2Dwave,[NaN;NaN;0]];

inputParameters = [inputParameters2Dwave,inputParametersRandom,[NaN;NaN;0]];

% controls(1:end-1) = controls2DWave;

controls(1:end-1) = [controls2DWave;controlsRandom];

controls{counter} = ones(m, N); % Constant input at the end


K = counter;

tStepEnd = toc(tStepStart);

fprintf(['done (',num2str(tStepEnd,3) ,'s) \n'])

fprintf(['There are ', num2str(K),' inputs in this experiment'])

%%

%%%% Collection of Fisher matrices

fprintf(['\n','Fisher Matrices computations... \n'])

% Use cell array for M
M = cell(K, 1);

tStepStart = tic;

for i = 1:K
    u = controls{i};
    u = u(:);
    Mi = zeros(p,p);
    L2u = max(1/Tf*(u'*u)/N, 0.000001);
    u = u/sqrt(L2u);
    
    Hju = zeros(N*q, p);
    SBiHju = zeros(N*q, p);

    for j = 1:p
        Hju(:,j) = Hb{j} * u;
        SBiHju(:,j) = SigmaBoldinv * Hju(:,j);
    end

    for j = 1:p
        for k = 1:j
            uLu_jk = (Hju(:,k)' * SBiHju(:,j) + Hju(:,j)' * SBiHju(:,k))/2;
            Mi(j,k) = uLu_jk + Q(j,k);
            Mi(k,j) = uLu_jk + Q(j,k);
        end
    end

    M{i} = Mi/(N*q); % Rescale of final Fisher infomatrix in order to account for influece of N and q
    
    if mod(i,floor(K/10))==0
        fprintf(['Progress:  %3.0f%% '],ceil(100*i/K))
        tStepEnd = toc(tStepStart);
        fprintf(['(',num2str(tStepEnd,3),' s).\n'])
        tStepStart = tic;
    end
end


tStepEnd = toc(tStepStart);

%% Input selection step

timeTakenCheck = toc(timeTakenStart);

fprintf(['\n','Total computation time up to now: ', num2str(timeTakenCheck,3) ,' s\n'])

fprintf([' \n'])

%%%% Finding the best vertex for initialization

maxScore = -Inf;
maxScoreIndex = 0;

for k = 1:K
    if maxScore < logdet(M{k})
        maxScore = logdet(M{k});
        maxScoreIndex = k;
    end
end

Mtop = M{maxScoreIndex};


maxBisec = 20;

wold = ones(K,1)/K;
wnew = wold;

% Compute weighted sum
Minit = zeros(p,p);
for k = 1:K
    Minit = Minit + wnew(k) * M{k};
end

Mold = Minit;
Mnew = Mold;

%%

maxStepSearch = 1000;
i = 0;
improvementThreshold = 10^(-8);
improvement = 1;

while (i < maxStepSearch) && (improvement > improvementThreshold)
    
    i = i+1;

    Mold = Mnew;
    wold = wnew;

    % Step 1: search for highest derivative

    maxDer = -Inf;
    maxDerIndex = 0;
    for k = 1:K                     
        der = trAB(inv(Mold), M{k});
        if maxDer < der
            maxDer = der;
            maxDerIndex = k;
        end
    end
    
    ks = maxDerIndex;
    
    % Step 2: search for best step

    Mks = M{ks};
    wks = zeros(K, 1);
    wks(ks) = 1;
    
    %%%% Adaptive step 
    
    % The adaptive step is computed via bisection in this method.
    
    a = 0;
    b = 1;
    
    % Assumption is that at a = 0, the function is positive and at b=1 it 
    % is negative. At a it makes sense by construction. However it may 
    % be positive at b in some cases. We don't care.
    
    for j = 1:maxBisec
        c = (a+b)/2;
        if trAB(inv((1-c)*Mold + c*Mks), (Mks-Mold)) >= 0
            a = c;
        else 
            b = c;
        end
    end
    alpha = (a+b)/2; % Final step
    
    
    Mnew = (1-alpha) * Mold + alpha * Mks;
    wnew = (1-alpha) * wold + alpha * wks;
    

    improvement = (logdet(Mnew)-logdet(Mold))/abs(logdet(Mold));
    
    if mod(i,10) == 0 
        fprintf(['Step of search : ', num2str(i), ', score improvement : ', num2str(improvement), '\n'])
    end

end



fprintf(['Max score of a matrix : ', num2str(maxScore), '\n'])

fprintf(['Final score : ', num2str(logdet(Mnew)), '\n'])

fprintf(['Score of averaged matrix : ', num2str(logdet(Minit)), '\n'])


threshold = 0.01;

idx = find(abs(wnew) > threshold);


timeTakenCheck = toc(timeTakenStart);

fprintf(['\n','Total computation time: ', num2str(timeTakenCheck,3) ,'s (',num2str(floor(timeTakenCheck/60)),' min ',num2str(mod(timeTakenCheck,60),2),'s)\n'])


%%

figure(1)
clf
pause(1)

% Define spatial domain
dx = retinalWidth/(sdn-1);
x = 0:dx:(sdn-1)*dx;

lidx = length(idx);
info = cell(lidx,4);
for k=1:length(idx)
    info{k,1} = num2str(wnew(idx(k)));
    info{k,2} = num2str(inputParameters(1,idx(k))/pi);
    info{k,3} = num2str(inputParameters(2,idx(k)));
    info{k,4} = num2str(inputParameters(3,idx(k)));
end

% Calculate subplot layout: 2 rows, variable columns
nCols = ceil(lidx/2);
nRows = min(2, lidx);  % Use 2 rows unless there's only 1 element

for l=1:N
    for k=1:length(idx)
        subplot(nRows, nCols, k)
        control_plot1D = controls{idx(k)}(:,l);
        control_plot = reshape(control_plot1D, sdn, sdn);
        
        imagesc(control_plot, [0,2])
        colormap(gray);
        axis image;
        
        xlabel('Position')
        ylabel('Position')
        title("weight = " + info{k,1} + ...
            ", \theta = " + info{k,2} + ...
            "\pi, c = " + info{k,3} + ...
            ", k = " + info{k,4})
        drawnow
        pause(0.001)
    end
end