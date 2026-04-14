% The following script reproduces the numerical experiments presented in:
%
%   "Optimal input design for parameter estimation in linear controlled
%    stochastic systems with retinal applications"
%   L. Sacchelli, B. Cessac, L. Pronzato
%   Preprint 2026
%
% GitHub repository: https://github.com/sacchelli/Retinal-experiment-design
%
% If you use this code in your work, please cite the above reference.
% This code is distributed under the MIT License (see LICENSE file).



addpath(genpath('./Functions'))
clc
timeTakenStart = tic;



%% FRONT MATTER %%


dn = 100; % Size of one layer$

N = 100; % Number of measurements


% Pick the rng seed if necessary.

seed = 1;

rng(seed);


%%% Should we cutoff correlations?

% If cutoff is 'true', we do not use the full matrices for correlations, 
% using the fact that there is a forgetting effect in the
% dynamics.

% cutoff = true;
cutoff = false;

% If cutoff is 'true', then depth determines the forgetting time
% beyond which we assume we can discard the correlations. It is
% possible to experiment with the function 'plotDiagonalMax' to see the
% decrease in amplitude of the matrices along the diagonal. Naturally, this
% effect is highly dependent on Delta. Taking Delta << tauX implies that
% cut off should be avoided, taking N Delta > 3 tau_X, it becomes natural
% to make this reduction.

depth = 20; 


%%% Should we add random controls in addition to waves?

addRandomInputs = true;

% How many?

nbrRandomInputs = 1000; % must be >0 if addRandomInputs is true.


%%% Should we create an animated gif of the result at the end?

%createAnimatedGIF = true;
createAnimatedGIF = false;


%% DYNAMICS CONSTRUCTION

%%%% Data // following the reference.

retinalWidth = 2; % mm 

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

n = 3 * dn; % size of the system

m = dn; % number of inputs

p = 7; % number of parameters

q = 4; % number of outputs

positions = [floor(dn/5),floor(2*dn/5),floor(3*dn/5),floor(4*dn/5)];


%positions = randi(dn,1,q);%[floor(dn/5),floor(2*dn/5),floor(3*dn/5),floor(4*dn/5)];



delta = retinalWidth/dn; % distance between cells


%%% Determination of size of the noise

% Notice that here, if we keep id as the noise matrix then only the
% relative sizes of the noises matter in the end.
% D is the dynamical noise and Sigmap is the output/measurement noise.

% D = .1 * eye(n); 
% Sigmap = 1 * eye(q); 

% D = .5 * eye(n); 
% Sigmap = 5 * eye(q); 

D = 1 * eye(n); 
Sigmap = 10 * eye(q); 

%gaussianPooling1DGraph(dn,delta,sigmaBG);

%%%% Model topology

% Below is the choice of topology for the B <-> A connection (matrix Gamma)
% All choices below are variations on the 1D Laplacian. Other choices
% could be made.

% 'adjacencyLaplacian' connection graph:
%           Bi
%         / | \
%        1  1  1
%       /   |   \
%   A(i-1)  Ai   A(i+1)

% 'closestNeighbor' connection graph:
%           Bi
%         / | \
%        1  0  1
%       /   |   \
%   A(i-1)  Ai   A(i+1)

% Two more possible examples

% 'absoluteLaplacian' connection graph:
%           Bi
%         / | \
%        1  2  1
%       /   |   \
%   A(i-1)  Ai   A(i+1)

% 'Laplacian' connection graph (unrealistic for the retina):
%           Bi
%         / | \
%        1  -2 1
%       /   |   \
%   A(i-1)  Ai   A(i+1)

Gamma = adjacencyLaplacian1DGraph(dn);  % Case 1-a

% Gamma = closestNeighbor1DGraph(dn); % Case 2-a

% Gamma = absoluteLaplacian1DGraph(dn);

% Gamma = laplacian1DGraph(dn);


% B -> G and A -> G connections are achieved by Gaussian pooling, (Gaussian
% kernel on the distance between nodes).

GammaBG = gaussianPooling1DGraph(dn,delta,sigmaBG);
GammaAG = gaussianPooling1DGraph(dn,delta,sigmaAG);

%%% Construction of the matrix A.

id = speye(dn);
zz = sparse(dn,dn);


Acell = cell(p+1, 1);

Acell{1} = [-id , zz, zz ; zz, zz, zz ; zz, zz, zz];
Acell{2} = [zz , zz, zz ; zz, -id, zz ; zz, zz, zz];
Acell{3} = [zz , zz, zz ; zz, zz, zz ; zz, zz, -id];

Acell{4} = [zz , - Gamma, zz ; zz, zz, zz ; zz, zz, zz];
Acell{5} = [zz , zz, zz ; Gamma, zz, zz ; zz, zz, zz];

Acell{6} = [zz , zz, zz ; zz, zz, zz ; GammaBG, zz, zz];
Acell{7} = [zz , zz, zz ; zz, zz, zz ; zz, - GammaAG, zz];

Acell{8} = 1/tauB*Acell{1} + 1/tauA*Acell{2} + 1/tauG*Acell{3} ...
    + wm*Acell{4} + wp*Acell{5} + wBG*Acell{6} + wAG*Acell{7};

% Acell{8} contains the experimental value, the others serve as variations
% around that value.

DeltaMin = 8.25; % ms

Delta = DeltaMin; % Time between two measurements

Tf = N * Delta; % Total duration of the experiment


%%% Input and output matrices

B = sparse([eye(dn,m);zeros(2*dn,m)]);

C0 = sparse([zeros(dn,2*dn),eye(dn,dn)]);

C = C0(positions,:);

%%% Contraction of A

thetaData = [1/tauB, 1/tauA, 1/tauG, wm, wp, wBG, wAG]; % Compilation of the experimental data

% We can choose the nominal vector theta, either as thetaData or as a
% variation of it. Again, we made the choice for this code that the 
% values input in thetaNominal are variations around thetaData (already 
% contained in Acell{8}).


thetaNominal = 0.*thetaData;

% thetaNominal = (rand(1,p)-1/2)/2.*thetaData;


A = Acell{p+1};
for i = 1:p
    A = A + thetaNominal(i) * Acell{i};
end

maxReEigA = max(real(eigs(A)));


fprintf(['Data created. Topology: Gamma(i,i)=',num2str(Gamma(1,1)),'. Size of one layer: ', num2str(dn),'. Length of time series: ', num2str(N),'. Max Re(λ(A)) = ', num2str(maxReEigA), '. \n'])


 
%% PRECOMPUTATION STEP %%

%%% Construction of sensitivities

% This construction follows the paper precisely.


fprintf(['\n','Dynamics and variance pre-computation... '])

tStepStart = tic;

[F,G] = buildFG(A, B, Delta);
Sigma = lyap(A,D*D');
SigmaDelta = Sigma-F*Sigma*F';

% Use a cutoff accelerate the computations by making the matrix banded:

if cutoff
    SigmaBoldInv = buildSigmaBoldInvCutoff(C, F, Sigma, Sigmap, N, depth);
else
    SigmaBoldInv = buildSigmaBoldInv(C, F, Sigma, Sigmap, N);
end


tStepEnd = toc(tStepStart);

fprintf(['done (', num2str(tStepEnd,3) ,'s). \n']) 


%% PARTIAL DERIVATIVES PRE-COMPUTATIONS %%

fprintf(['\n','Parameter sensitivities pre-computation... \n'])

Hb = cell(p, 1);
dSb = cell(p, 1);

for i = 1:p
    fprintf(['Parameter ',num2str(i),'/',num2str(p),' '])
    tStepStart = tic;
    
    Ab = [A, zeros(n,n); Acell{i}, A];
    
    Bb = [B; zeros(n,m)];
    Cb = [zeros(q,n), C];
    Cbp = [C, zeros(q,n)];
    Db = [D, zeros(n,n); zeros(n,n), zeros(n,n)];

    [Fb,Gb] = buildFG(Ab, Bb, Delta);
    Sigmab = lyap(Ab,Db*Db'); 
    
    if cutoff 
        dSb{i} = buildSigmaBoldCutoff(Cb, Cbp, Fb, Sigmab, zeros(q,q), N, depth);
        Hb{i} = buildHCutoff(Cb, Fb, Gb, N, depth);
    else
        dSb{i} = buildSigmaBold(Cb, Cbp, Fb, Sigmab, zeros(q,q), N);
        Hb{i} = buildH(Cb, Fb, Gb, N);
    end 
    
    tStepEnd = toc(tStepStart);

    fprintf(['done (', num2str(tStepEnd,3) ,'s). \n'])
end


%% COMPUTATION OF MATRIX Q %%

timeTakenCheck = toc(timeTakenStart);

fprintf(['\n','Total computation time up to now: ', num2str(timeTakenCheck,3) ,'s\n'])

fprintf(['\n','Q matrix computation... \n'])

if cutoff
    Q = buildQCutoff(SigmaBoldInv, dSb, N, depth);
else
    Q = buildQ(SigmaBoldInv, dSb, N);
end

timeTakenCheck = toc(timeTakenStart);

fprintf(['\n','Q matrix computation done. \n'])

fprintf(['\n','Total computation time up to now: ', num2str(timeTakenCheck,3) ,'s\n'])


%% COLLECTION OF INPUTS %%

%%% Test family: waves 

fprintf(['\n','Input collection... '])

tStepStart = tic;

% cStar = retinalWidth/1000;

cStar = retinalWidth/(100*Delta);
% cStar = retinalWidth/(Tf);

cMax = 25*cStar;   % max wave velocity
cMin = 0; % min wave velocity

kMax = pi/retinalWidth*dn/2;   % max wave number
kMin = pi/retinalWidth/10;   % min wave number

cN = 101;
kN = 100;

[controls1DWave, inputParameters1Dwave] = generate1DPlaneWaveControls(dn, N, Delta, Tf, retinalWidth, cMin, cMax, cN, kMin, kMax, kN);

if addRandomInputs % produce random inputs if asked
    %[controlsRandom, inputParametersRandom] = generate1DRandomControls(dn, N, nbrRandomInputs, seed);
    [controlsRandom, inputParametersRandom] = generate1DRandomBinaryControls(dn, N, nbrRandomInputs, seed);
end


% Concatenate inputs and constant input.
if addRandomInputs 
    counter = length(controls1DWave) + length(controlsRandom) + 1;
    inputParameters = [inputParameters1Dwave, inputParametersRandom, [NaN;0]];
    controls = cell(counter,1);
    controls(1:end-1) = [controls1DWave; controlsRandom];
else
    counter = length(controls1DWave) + 1;
    inputParameters = [inputParameters1Dwave,  [NaN;0]];
    controls = cell(counter,1);
    controls(1:end-1) = controls1DWave;
end

%

controls{counter} = ones(m, N); % Constant input at the end


K = counter;


for i = 1:K
    uCol = controls{i}(:);
    L2u = max((uCol'*uCol)/(N), 0.000001);
    controls{i} = controls{i}/sqrt(L2u);
end



tStepEnd = toc(tStepStart);

fprintf(['done (',num2str(tStepEnd,3) ,'s) \n'])

fprintf(['There are ', num2str(K),' inputs in this experiment'])



%% FISHER MATRICES CONSTRUCTION %%

% Please notice that the Fisher matrices are rescaled 

%%% Collection of Fisher matrices

fprintf(['\n','Fisher Matrices computations... \n'])

M = cell(K, 1);

tStepStart = tic;

Mi = zeros(p,p);
Hju = zeros(N*q, p);
SBiHju = zeros(N*q, p);


timebar = waitbar(0, 'FIM computation');
tLoopStart = tic;

for i = 1:K

    if mod(i,50)==0
        timeInLoop = toc(tLoopStart);
        timeRemainingMin=floor((K-i)/i*timeInLoop/60);
        timeRemainingSec=floor((K-i)/i*timeInLoop-60*timeRemainingMin);
        
        waitbar(i/K, timebar, ['FIM computation (',...
        num2str(floor(i/K*100)),' %, est. remaining time: ',...
        num2str(timeRemainingMin),'min ',...
        num2str(timeRemainingSec),'s)']);
    end


    u = controls{i};
    u = u(:);
    %L2u = max((u'*u)/N, 0.000001);
    %u = u/sqrt(L2u);
    
    for j = 1:p
        Hju(:,j) = Hb{j} * u;
        SBiHju(:,j) = SigmaBoldInv * Hju(:,j);
    end

    for j = 1:p
        for k = 1:j
            uLu_jk = (Hju(:,k)' * SBiHju(:,j) + Hju(:,j)' * SBiHju(:,k))/2;
            Mi(j,k) = uLu_jk + Q(j,k);
            Mi(k,j) = uLu_jk + Q(j,k);
        end
    end

    M{i} = Mi; 
    
    if mod(i,floor(K/10))==0
        fprintf(['Progress:  %3.0f%% '],ceil(100*i/K))
        tStepEnd = toc(tStepStart);
        fprintf(['(',num2str(tStepEnd,3),' s).\n'])
        tStepStart = tic;
    end
end

close(timebar);

tStepEnd = toc(tStepStart);

%% D-OPTIMALITY %%

timeTakenCheck = toc(timeTakenStart);

fprintf(['\n','Total computation time up to now: ', num2str(timeTakenCheck,3) ,' s\n'])

fprintf([' \n'])

%%% Finding the best vertex before initialization

maxScore = -Inf;
maxScoreIndex = 0;

for k = 1:K
    if maxScore < logdet(M{k})
        maxScore = logdet(M{k});
        maxScoreIndex = k;
    end
end

minScore = Inf;
minScoreIndex = 0;

for k = 1:K
    if minScore > logdet(M{k})
        minScore = logdet(M{k});
        minScoreIndex = k;
    end
end


if maxScoreIndex>0
    Mtop = M{maxScoreIndex};
else
    Mtop = NaN;
end

maxBisec = 20;  % How far should we push bisection.

%%% Initialize at the average input

wold = ones(K,1)/K;
wnew = wold;



% Compute weighted sum
Minit = zeros(p,p);
for k = 1:K
    Minit = Minit + wnew(k) * M{k};
end

Mold = Minit;
Mnew = Mold;

%%% Exploration loop

maxStepSearch = 10000; % Max number of loops
i = 0;
improvementThreshold = 10^(-8);
improvement = 1;

maxDerThreshold = 1;
maxDer = Inf;

tSearchStart = tic;

while (i < maxStepSearch) && (improvement > improvementThreshold)
    
    i = i+1;

    Mold = Mnew;
    wold = wnew;

    % Step 1: search for highest derivative

    maxDer = -Inf;
    maxDerIndex = 0;

    
    
    MoldInv = Mold\eye(p);

    for k = 1:K                     
        der = trAB(MoldInv, M{k});
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

    if mod(i,100) == 0 
        fprintf(['Step of search : ', num2str(i), ', score improvement: ', num2str(improvement), '\n'])
    end

    
    % if mod(i,100) == 0 
    %     fprintf(['Step of search : ', num2str(i), ', max derivative: ', num2str(maxDer), ', index: ',num2str(maxDerIndex),'\n'])
    % end
end


MoptD = Mnew;
wnew_D = wnew;

tSearchTaken = toc(tSearchStart);

fprintf(['\n','Search time: ', num2str(tSearchTaken,3) ,'s \n'])

fprintf(['Max score of a matrix : ', num2str(maxScore), '\n'])

fprintf(['Final score : ', num2str(logdet(Mnew)), '\n'])

fprintf(['Score of averaged matrix : ', num2str(logdet(Minit)), '\n'])


if addRandomInputs
    wRand = 1/nbrRandomInputs;
    Mrand = zeros(p,p);
    for k = K-1:-1:K-nbrRandomInputs
        Mrand = Mrand + wRand * M{k};
    end
    fprintf(['Score of random inputs : ', num2str(logdet(Mrand)), '\n'])
end


threshold = 0.01;

idx = find(abs(wnew) > threshold);

timeTakenCheck = toc(timeTakenStart);

fprintf(['\n','Total computation time: ', num2str(timeTakenCheck,3) ,'s (',num2str(floor(timeTakenCheck/60)),' min ',num2str(mod(timeTakenCheck,60),2),'s)\n'])





%% A-OPTIMALITY %%

timeTakenCheck = toc(timeTakenStart);

fprintf(['\n','Total computation time up to now: ', num2str(timeTakenCheck,3) ,' s\n'])

fprintf([' \n'])

%%% Finding the best vertex before initialization

maxScore = -Inf;
maxScoreIndex = 0;

Phi = @(M) 1/trace(inv(M));


for k = 1:K
    if maxScore < Phi(M{k})
        maxScore = Phi(M{k});
        maxScoreIndex = k;
    end
end


if maxScoreIndex>0
    Mtop = M{maxScoreIndex};
else
    Mtop = NaN;
end

maxBisec = 20;  % How far should we push bisection.

%%% Initialize at the average input

wold = ones(K,1)/K;
wnew = wold;



% Compute weighted sum
Minit = zeros(p,p);
for k = 1:K
    Minit = Minit + wnew(k) * M{k};
end

Mold = Minit;
Mnew = Mold;

%%% Exploration loop

maxStepSearch = 10000; % Max number of loops
i = 0;
improvementThreshold = 10^(-8);
improvement = 1;

maxDerThreshold = 1;
maxDer = Inf;

tSearchStart = tic;

while (i < maxStepSearch) && (improvement > improvementThreshold)
    
    i = i+1;

    Mold = Mnew;
    wold = wnew;

    % Step 1: search for highest derivative

    maxDer = -Inf;
    maxDerIndex = 0;

    
    
    MoldInv = Mold\eye(p);
    PhiM = Phi(Mold);
    for k = 1:K                     
        der =  PhiM^2 * trAB(MoldInv*MoldInv, M{k});
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
        Mc = (1-c)*Mold + c*Mks;
        McInv = Mc\eye(p);
        gradPhiMc = Phi(Mc)^2 * McInv*McInv;

        if trAB(gradPhiMc, (Mks-Mold)) >= 0
            a = c;
        else 
            b = c;
        end
    end
    alpha = (a+b)/2; % Final step
    
    
    Mnew = (1-alpha) * Mold + alpha * Mks;
    wnew = (1-alpha) * wold + alpha * wks;
    

    improvement = (Phi(Mnew)-Phi(Mold))/abs(Phi(Mold));

    if mod(i,100) == 0 
        fprintf(['Step of search : ', num2str(i), ', score improvement: ', num2str(improvement), '\n'])
    end

    
    % if mod(i,100) == 0 
    %     fprintf(['Step of search : ', num2str(i), ', max derivative: ', num2str(maxDer), ', index: ',num2str(maxDerIndex),'\n'])
    % end
end


MoptA = Mnew;
wnew_A = wnew;




tSearchTaken = toc(tSearchStart);

fprintf(['\n','Search time: ', num2str(tSearchTaken,3) ,'s \n'])

fprintf(['Max score of a matrix : ', num2str(maxScore), '\n'])

fprintf(['Final score : ', num2str(Phi(Mnew)), '\n'])

fprintf(['Score of averaged matrix : ', num2str(Phi(Minit)), '\n'])


if addRandomInputs
    wRand = 1/nbrRandomInputs;
    Mrand = zeros(p,p);
    for k = K-1:-1:K-nbrRandomInputs
        Mrand = Mrand + wRand * M{k};
    end
    fprintf(['Score of random inputs : ', num2str(Phi(Mrand)), '\n'])
end


threshold = 0.01;

idx = find(abs(wnew) > threshold);

timeTakenCheck = toc(timeTakenStart);

fprintf(['\n','Total computation time: ', num2str(timeTakenCheck,3) ,'s (',num2str(floor(timeTakenCheck/60)),' min ',num2str(mod(timeTakenCheck,60),2),'s)\n'])





%% C-OPTIMALITY %%

timeTakenCheck = toc(timeTakenStart);

fprintf(['\n','Total computation time up to now: ', num2str(timeTakenCheck,3) ,' s\n'])

fprintf([' \n'])

%%% Finding the best vertex before initialization

maxScore = -Inf;
maxScoreIndex = 0;

P = zeros(p,1);
P(3) = 1;

Phic = @(M) 1/trace(P'*(M\P));


for k = 1:K
    if maxScore < Phic(M{k})
        maxScore = Phic(M{k});
        maxScoreIndex = k;
    end
end


if maxScoreIndex>0
    Mtop = M{maxScoreIndex};
else
    Mtop = NaN;
end

maxBisec = 20;  % How far should we push bisection.

%%% Initialize at the average input

wold = ones(K,1)/K;
wnew = wold;



% Compute weighted sum
Minit = zeros(p,p);
for k = 1:K
    Minit = Minit + wnew(k) * M{k};
end

Mold = Minit;
Mnew = Mold;

%%% Exploration loop

maxStepSearch = 10000; % Max number of loops
i = 0;
improvementThreshold = 10^(-8);
improvement = 1;

maxDerThreshold = 1;
maxDer = Inf;

tSearchStart = tic;

while (i < maxStepSearch) && (improvement > improvementThreshold)
    
    i = i+1;

    Mold = Mnew;
    wold = wnew;

    % Step 1: search for highest derivative

    maxDer = -Inf;
    maxDerIndex = 0;

    
    
    MoldInv = Mold\eye(p);
    PhiM = Phic(Mold);
    for k = 1:K                     
        der =  PhiM^2 * trAB(MoldInv*(P*P')*MoldInv, M{k});
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
        Mc = (1-c)*Mold + c*Mks;
        McInv = Mc\eye(p);
        gradPhiMc = Phic(Mc)^2 * McInv*(P*P')*McInv;

        if trAB(gradPhiMc, (Mks-Mold)) >= 0
            a = c;
        else 
            b = c;
        end
    end
    alpha = (a+b)/2; % Final step
    
    
    Mnew = (1-alpha) * Mold + alpha * Mks;
    wnew = (1-alpha) * wold + alpha * wks;
    

    improvement = (Phic(Mnew)-Phic(Mold))/abs(Phic(Mold));

    if mod(i,100) == 0 
        fprintf(['Step of search : ', num2str(i), ', score improvement: ', num2str(improvement), '\n'])
    end

    
    % if mod(i,100) == 0 
    %     fprintf(['Step of search : ', num2str(i), ', max derivative: ', num2str(maxDer), ', index: ',num2str(maxDerIndex),'\n'])
    % end
end


MoptC = Mnew;
wnew_C = wnew;




tSearchTaken = toc(tSearchStart);

fprintf(['\n','Search time: ', num2str(tSearchTaken,3) ,'s \n'])

fprintf(['Max score of a matrix : ', num2str(maxScore), '\n'])

fprintf(['Final score : ', num2str(Phic(Mnew)), '\n'])

fprintf(['Score of averaged matrix : ', num2str(Phic(Minit)), '\n'])


if addRandomInputs
    wRand = 1/nbrRandomInputs;
    Mrand = zeros(p,p);
    for k = K-1:-1:K-nbrRandomInputs
        Mrand = Mrand + wRand * M{k};
    end
    fprintf(['Score of random inputs : ', num2str(Phic(Mrand)), '\n'])
end

Mrand_nominal = Mrand

threshold = 0.01;

idx = find(abs(wnew) > threshold);

timeTakenCheck = toc(timeTakenStart);

fprintf(['\n','Total computation time: ', num2str(timeTakenCheck,3) ,'s (',num2str(floor(timeTakenCheck/60)),' min ',num2str(mod(timeTakenCheck,60),2),'s)\n'])



%% representation of theoretical gains

Nexp = 100;

figure(1)
clf


% Define parameter names
param_names = {'\tau_B^{-1}', '\tau_A^{-1}', '\tau_G^{-1}', 'w^-', 'w^+', 'w^B_G','w^A_G'};
param_names_plain = {'1/tau_B', '1/tau_A', '1/tau_G', 'w^-', 'w^+', 'w^B_G','w^A_G'};


% --- Subplot 1: Eigenvalues (Comparison) ---
subplot(1,3,1)
% Use semilogy to match the style of your other eigenvalue plot
semilogy(1:p,1/sqrt(Nexp)*sqrt(eig(inv(MoptD))), 'b-o', 'LineWidth', 2, 'MarkerSize', 8)
hold on
%semilogy(1:p, 1/sqrt(Nexp)*sqrt(eig(inv(Mrand))), 'r-s', 'LineWidth', 2, 'MarkerSize', 8)
semilogy(1:p,1/sqrt(Nexp)*sqrt(eig(inv(MoptA))), 'g-o', 'LineWidth', 2, 'MarkerSize', 8)
semilogy(1:p,1/sqrt(Nexp)*sqrt(eig(inv(MoptC))), 'm-o', 'LineWidth', 2, 'MarkerSize', 8)
grid on
xlim([1,p])
title('Predicted Uncertainty ellipsoid')
ylabel('\surd eigenvalues')
xlabel('Eigenvalue Index')
legend('D-optimal', ...%'Random', ...
    'A-optimal','c-optimal','Location', 'best')
set(gca, 'FontSize', 11)



% --- Subplot 2: Parameter uncertainties (square root of diagonal of Inverse) ---
subplot(1,3,2)

% Plotting diag(inv(M)) gives the variance of each parameter
plot(1:p, 1/sqrt(Nexp)*sqrt(diag(inv(MoptD))), 'b-o', 'LineWidth', 2, 'MarkerSize', 8)
hold on
%plot(1:p, 1/sqrt(Nexp)*sqrt(diag(inv(Mrand))), 'r-s', 'LineWidth', 2, 'MarkerSize', 8)
plot(1:p, 1/sqrt(Nexp)*sqrt(diag(inv(MoptA))), 'g-o', 'LineWidth', 2, 'MarkerSize', 8)
plot(1:p, 1/sqrt(Nexp)*sqrt(diag(inv(MoptC))), 'm-o', 'LineWidth', 2, 'MarkerSize', 8)
%plot(1:p, 1/sqrt(3)*1/sqrt(Nexp)*sqrt(diag(inv(Mrand))), 'y-.s', 'LineWidth', 2, 'MarkerSize', 8)
grid on
xlim([1,p])
ymax = max([max(1/sqrt(Nexp)*sqrt(diag(inv(MoptD)))), max(1/sqrt(Nexp)*sqrt(diag(inv(MoptA)))), max(1/sqrt(Nexp)*sqrt(diag(inv(MoptC))))]);

ymax = ceil(20*10^(-ceil(log10(ymax)))*ymax)/(20*10^(-ceil(log10(ymax))));



ylim([0,ymax])
title('Predicted parameter uncertainty')
ylabel('StdDev(\theta_i)')
xlabel('Parameter Index')
legend('D-optimal', ...%'Random',...
     'A-optimal','c-optimal','Location', 'best')
set(gca, 'FontSize', 11)
set(gca, 'XTick', 1:p, 'XTickLabel', param_names, 'TickLabelInterpreter', 'tex')






% --- Subplot 3: Parameter Variances (Diagonal of Inverse) ---
subplot(1,3,3)

% Plotting diag(inv(M)) gives the variance of each parameter
plot(1:p, 1/sqrt(Nexp)*sqrt(diag(inv(MoptD)))./thetaData', 'b-o', 'LineWidth', 2, 'MarkerSize', 8)
hold on
%plot(1:p, 1/sqrt(Nexp)*sqrt(diag(inv(Mrand)))./thetaData', 'r-s', 'LineWidth', 2, 'MarkerSize', 8)
plot(1:p, 1/sqrt(Nexp)*sqrt(diag(inv(MoptA)))./thetaData', 'g-o', 'LineWidth', 2, 'MarkerSize', 8)
plot(1:p, 1/sqrt(Nexp)*sqrt(diag(inv(MoptC)))./thetaData', 'm-o', 'LineWidth', 2, 'MarkerSize', 8)
%plot(1:p, 1/sqrt(3)*1/sqrt(Nexp)*sqrt(diag(inv(Mrand)))./thetaData', 'y-.s', 'LineWidth', 2, 'MarkerSize', 8)
grid on
xlim([1,p])

ymax = max([max(1/sqrt(Nexp)*sqrt(diag(inv(MoptD)))./thetaData'),max(1/sqrt(Nexp)*sqrt(diag(inv(MoptA)))./thetaData'),max(1/sqrt(Nexp)*sqrt(diag(inv(MoptC)))./thetaData')]);

ymax = ceil(20*ymax)/20;

ylim([0,ymax])
title('Relative predicted parameter uncertainty')
ylabel('StdDev(\theta_i)')
xlabel('Parameter Index')
legend('D-optimal', ...%'Random',...
     'A-optimal','c-optimal','Location', 'best')
set(gca, 'FontSize', 11)
set(gca, 'XTick', 1:p, 'XTickLabel', param_names, 'TickLabelInterpreter', 'tex')


sgtitle(sprintf('Predicted uncertainty for D-Optimal vs Random Inputs  (K=%d)', Nexp), ...
        'FontSize', 16, 'FontWeight', 'bold')




%% representation of theoretical gains

Nexp = 100;

figure(2)
clf


% Define parameter names
param_names = {'\tau_B^{-1}', '\tau_A^{-1}', '\tau_G^{-1}', 'w^-', 'w^+', 'w^B_G','w^A_G'};
param_names_plain = {'1/tau_B', '1/tau_A', '1/tau_G', 'w^-', 'w^+', 'w^B_G','w^A_G'};


% --- Subplot 1: Eigenvalues (Comparison) ---
subplot(1,2,1)
% Use semilogy to match the style of your other eigenvalue plot
semilogy(1:p,1/sqrt(Nexp)*sqrt(eig(inv(MoptD))), 'b-o', 'LineWidth', 2, 'MarkerSize', 8)
hold on
%semilogy(1:p, 1/sqrt(Nexp)*sqrt(eig(inv(Mrand))), 'r-s', 'LineWidth', 2, 'MarkerSize', 8)
semilogy(1:p,1/sqrt(Nexp)*sqrt(eig(inv(MoptA))), 'g-o', 'LineWidth', 2, 'MarkerSize', 8)
semilogy(1:p,1/sqrt(Nexp)*sqrt(eig(inv(MoptC))), 'm-o', 'LineWidth', 2, 'MarkerSize', 8)
grid on
xlim([1,p])
title('Predicted Uncertainty ellipsoid')
ylabel('\surd eigenvalues')
xlabel('Eigenvalue Index')
legend('D-optimal', ...%'Random', ...
    'A-optimal','c-optimal','Location', 'best')
set(gca, 'FontSize', 11)



% --- Subplot 2: Parameter Variances (Diagonal of Inverse) ---
subplot(1,2,2)

% Plotting diag(inv(M)) gives the variance of each parameter
plot(1:p, 1/sqrt(Nexp)*sqrt(diag(inv(MoptD))), 'b-o', 'LineWidth', 2, 'MarkerSize', 8)
hold on
%plot(1:p, 1/sqrt(Nexp)*sqrt(diag(inv(Mrand))), 'r-s', 'LineWidth', 2, 'MarkerSize', 8)
plot(1:p, 1/sqrt(Nexp)*sqrt(diag(inv(MoptA))), 'g-o', 'LineWidth', 2, 'MarkerSize', 8)
plot(1:p, 1/sqrt(Nexp)*sqrt(diag(inv(MoptC))), 'm-o', 'LineWidth', 2, 'MarkerSize', 8)
%plot(1:p, 1/sqrt(3)*1/sqrt(Nexp)*sqrt(diag(inv(Mrand))), 'y-.s', 'LineWidth', 2, 'MarkerSize', 8)
grid on
xlim([1,p])
ymax = max([max(1/sqrt(Nexp)*sqrt(diag(inv(MoptD)))), max(1/sqrt(Nexp)*sqrt(diag(inv(MoptA)))), max(1/sqrt(Nexp)*sqrt(diag(inv(MoptC))))]);

ymax = ceil(20*10^(-ceil(log10(ymax)))*ymax)/(20*10^(-ceil(log10(ymax))));



ylim([0,ymax])
title('Predicted parameter uncertainty')
ylabel('StdDev(\theta_i)')
xlabel('Parameter Index')
legend('D-optimal', ...%'Random',...
     'A-optimal','c-optimal','Location', 'best')
set(gca, 'FontSize', 11)
set(gca, 'XTick', 1:p, 'XTickLabel', param_names, 'TickLabelInterpreter', 'tex')




%% ROBUSTNESS CHECK: RE-EVALUATE OPTIMAL DESIGNS AT A NEIGHBOURING THETA %%
%
% The optimal designs (weights wnew_D, wnew_A, wnew_C) were found under the
% nominal theta. Here we ask: if the true parameter is actually theta_neighbor,
% what Fisher matrix do those same optimal inputs actually deliver?
%
% No new optimisation is needed. We only need to:
%   1. Build A_nbr from theta_neighbor and recompute F, G, Sigma, SigmaBoldInv.
%   2. Recompute H{i} and Q for theta_neighbor (same Acell{i} structure).
%   3. Recompute M_nbr{k} for every input k using the new quantities.
%   4. Form the weighted sums with wnew_D / wnew_A / wnew_C and compare.
%
% PREREQUISITES -- add these lines at the end of each optimisation section:
%   wnew_D = wnew;          % after D-optimal loop
%   wnew_A = wnew;          % after A-optimal loop
%   wnew_C = wnew;          % after C-optimal loop
%   Mrand_nominal = Mrand;  % right after the first Mrand computation (D section)

% -------------------------------------------------------------------------
% 1.  Define theta_neighbor
% -------------------------------------------------------------------------

perturbationFactor = 0.5;   % 20% off -- change sign or magnitude freely

% Multiplicative perturbation on every parameter simultaneously.
% To perturb individual entries instead, e.g.:
%   thetaNeighbor = thetaData;  thetaNeighbor(3) = thetaData(3) * 1.20;


rng(7);

thetaNeighbor    = thetaData .* (1 + perturbationFactor);

thetaNeighbor    = thetaData .* (1 + 2*perturbationFactor*(rand(1,p)-.5));

thetaNominal_nbr = thetaNeighbor - thetaData;   % additive shift on top of Acell{8}

fprintf(['\n', repmat('-',1,60), '\n'])
fprintf('ROBUSTNESS CHECK  --  theta_neighbor = theta_true * (1 + %.0f%%)\n', ...
        perturbationFactor*100)
fprintf([repmat('-',1,60), '\n'])

% -------------------------------------------------------------------------
% 2.  Build A_nbr and precompute dynamics quantities
% -------------------------------------------------------------------------

fprintf('\nDynamics pre-computation at theta_neighbor... ')
tStepStart = tic;

A_nbr = Acell{p+1};
for i = 1:p
    A_nbr = A_nbr + thetaNominal_nbr(i) * Acell{i};
end

B_nbr=B;

[F_nbr, ~] = buildFG(A_nbr, B_nbr, Delta);
Sigma_nbr = lyap(A_nbr, D*D');

if cutoff
    SigmaBoldInv_nbr = buildSigmaBoldInvCutoff(C, F_nbr, Sigma_nbr, Sigmap, N, depth);
else
    SigmaBoldInv_nbr = buildSigmaBoldInv(C, F_nbr, Sigma_nbr, Sigmap, N);
end

fprintf(['done (', num2str(toc(tStepStart), 3), 's).\n'])

% -------------------------------------------------------------------------
% 3.  Recompute H and Q at theta_neighbor
% -------------------------------------------------------------------------

fprintf('\nParameter sensitivities at theta_neighbor...\n')

Hb_nbr  = cell(p, 1);
dSb_nbr = cell(p, 1);

for i = 1:p
    fprintf('  Parameter %d/%d ', i, p)
    tStepStart = tic;

    Ab_nbr  = [A_nbr,     zeros(n,n);
               Acell{i},  A_nbr     ];
    Bb_nbr  = [B;          zeros(n,m)];
    Cb_nbr  = [zeros(q,n), C        ];
    Cbp_nbr = [C,          zeros(q,n)];
    Db_nbr  = [D,          zeros(n,n);
               zeros(n,n), zeros(n,n)];

    [Fb_nbr, Gb_nbr] = buildFG(Ab_nbr, Bb_nbr, Delta);
    Sigmab_nbr  = lyap(Ab_nbr, Db_nbr * Db_nbr');

    if cutoff
        dSb_nbr{i} = buildSigmaBoldCutoff(Cb_nbr, Cbp_nbr, Fb_nbr, Sigmab_nbr, zeros(q,q), N, depth);
        Hb_nbr{i}  = buildHCutoff(Cb_nbr, Fb_nbr, Gb_nbr, N, depth);
    else
        dSb_nbr{i} = buildSigmaBold(Cb_nbr, Cbp_nbr, Fb_nbr, Sigmab_nbr, zeros(q,q), N);
        Hb_nbr{i}  = buildH(Cb_nbr, Fb_nbr, Gb_nbr, N);
    end

    fprintf(['done (', num2str(toc(tStepStart), 3), 's).\n'])
end

fprintf('\nQ matrix at theta_neighbor... ')
tStepStart = tic;

if cutoff
    Q_nbr = buildQCutoff(SigmaBoldInv_nbr, dSb_nbr, N, depth);
else
    Q_nbr = buildQ(SigmaBoldInv_nbr, dSb_nbr, N);
end

fprintf(['done (', num2str(toc(tStepStart), 3), 's).\n'])

% -------------------------------------------------------------------------
% 4.  Re-evaluate Fisher matrices at theta_neighbor for all inputs
% -------------------------------------------------------------------------

fprintf('\nFisher matrices at theta_neighbor...\n')
tStepStart = tic;

M_nbr      = cell(K, 1);
Mi_nbr     = zeros(p, p);
Hju_nbr    = zeros(N*q, p);
SBiHju_nbr = zeros(N*q, p);


comp=zeros(K,1);
for k = K-1:-1:K-nbrRandomInputs
    comp(k)=1;
end

comp=comp+((abs(wnew_D)+abs(wnew_A)+abs(wnew_C))>0.000000001);


for i = 1:K
    u = controls{i}(:);
    %if comp(i)>0
        for j = 1:p
            Hju_nbr(:,j)    = Hb_nbr{j} * u;
            SBiHju_nbr(:,j) = SigmaBoldInv_nbr * Hju_nbr(:,j);
        end
    
        for j = 1:p
            for k = 1:j
                uLu_jk      = (Hju_nbr(:,k)' * SBiHju_nbr(:,j) + ...
                               Hju_nbr(:,j)' * SBiHju_nbr(:,k)) / 2;
                Mi_nbr(j,k) = uLu_jk + Q_nbr(j,k);
                Mi_nbr(k,j) = Mi_nbr(j,k);
            end
        end
    
        M_nbr{i} = Mi_nbr;
    %else
    %    M_nbr{i} = zeros(p);
    %end
    if mod(i,floor(K/10))==0
        fprintf(['Progress:  %3.0f%% '],ceil(100*i/K))
        tStepEnd = toc(tStepStart);
        fprintf(['(',num2str(tStepEnd,3),' s).\n'])
        tStepStart = tic;
    end
end

fprintf(['done (', num2str(toc(tStepStart), 3), 's).\n'])

% -------------------------------------------------------------------------
% 5.  Form weighted FIMs at theta_neighbor using the NOMINAL weights,
%     and re-evaluate the uniform random baseline at theta_neighbor
% -------------------------------------------------------------------------

MoptD_nbr = zeros(p,p);
MoptA_nbr = zeros(p,p);
MoptC_nbr = zeros(p,p);
for k = 1:K
    MoptD_nbr = MoptD_nbr + wnew_D(k) * M_nbr{k};
    MoptA_nbr = MoptA_nbr + wnew_A(k) * M_nbr{k};
    MoptC_nbr = MoptC_nbr + wnew_C(k) * M_nbr{k};
end

% Uniform random baseline re-evaluated at theta_neighbor
wRand     = 1 / nbrRandomInputs;
Mrand_nbr = zeros(p,p);
for k = K-1:-1:K-nbrRandomInputs
    Mrand_nbr = Mrand_nbr + wRand * M_nbr{k};
end

fprintf('\n--- Scores: nominal vs neighbour theta ---\n')
fprintf('                    nominal        neighbour\n')
fprintf('D-design  logdet:   %10.4f     %10.4f\n', logdet(MoptD),         logdet(MoptD_nbr))
fprintf('A-design  Phi_A:    %10.4f     %10.4f\n', Phi(MoptA),            Phi(MoptA_nbr))
fprintf('C-design  Phi_C:    %10.4f     %10.4f\n', Phic(MoptC),           Phic(MoptC_nbr))
fprintf('Random    logdet:   %10.4f     %10.4f\n', logdet(Mrand_nominal), logdet(Mrand_nbr))



% ------------------------------------------------------
%5'. 
% D-OPTIMALITY AT THETA_NEIGHBOR %%
% We already have M_nbr{k} for all k, so we just re-run the same
% Frank-Wolfe loop as the nominal D-optimal section, but on M_nbr.
% This gives the TRUE best achievable performance at theta_neighbor,
% against which we can compare MoptD_nbr (the nominal design evaluated
% at the neighbour).

fprintf(['\n', repmat('-',1,60), '\n'])
fprintf('D-optimal search at theta_neighbor...\n')
fprintf([repmat('-',1,60), '\n'])

%%% Find best vertex as starting upper bound
maxScore_nbr = -Inf;
for k = 1:K
    if ~isempty(M_nbr{k}) && logdet(M_nbr{k}) > maxScore_nbr
        maxScore_nbr = logdet(M_nbr{k});
    end
end

%%% Initialise at uniform weights
wold_nbr = ones(K,1)/K;
wnew_nbr = wold_nbr;

Mnew_nbr = zeros(p,p);
for k = 1:K
    Mnew_nbr = Mnew_nbr + wnew_nbr(k) * M_nbr{k};
end
Mold_nbr = Mnew_nbr;

%%% Frank-Wolfe loop
maxStepSearch    = 10000;
improvementThreshold = 1e-8;
improvement_nbr  = 1;
i_nbr            = 0;
maxBisec         = 20;

tSearchStart_nbr = tic;

while (i_nbr < maxStepSearch) && (improvement_nbr > improvementThreshold)

    i_nbr    = i_nbr + 1;
    Mold_nbr = Mnew_nbr;
    wold_nbr = wnew_nbr;

    % Step 1: find input with highest directional derivative
    MoldInv_nbr  = Mold_nbr \ eye(p);
    maxDer_nbr   = -Inf;
    maxDerIdx_nbr = 0;

    for k = 1:K
        if ~isempty(M_nbr{k})
            der = trAB(MoldInv_nbr, M_nbr{k});
            if der > maxDer_nbr
                maxDer_nbr    = der;
                maxDerIdx_nbr = k;
            end
        end
    end

    % Step 2: bisection line search
    Mks_nbr  = M_nbr{maxDerIdx_nbr};
    wks_nbr  = zeros(K,1);
    wks_nbr(maxDerIdx_nbr) = 1;

    a = 0;  b = 1;
    for j = 1:maxBisec
        c  = (a+b)/2;
        Mc = (1-c)*Mold_nbr + c*Mks_nbr;
        if trAB(Mc\eye(p), (Mks_nbr - Mold_nbr)) >= 0
            a = c;
        else
            b = c;
        end
    end
    alpha_nbr = (a+b)/2;

    Mnew_nbr = (1-alpha_nbr)*Mold_nbr + alpha_nbr*Mks_nbr;
    wnew_nbr = (1-alpha_nbr)*wold_nbr + alpha_nbr*wks_nbr;

    improvement_nbr = (logdet(Mnew_nbr) - logdet(Mold_nbr)) / abs(logdet(Mold_nbr));

    if mod(i_nbr, 100) == 0
        fprintf('Step %d, improvement: %e\n', i_nbr, improvement_nbr)
    end
end

MoptD_nbr_true = Mnew_nbr;   % best achievable at theta_neighbor under D-criterion
wnew_D_nbr     = wnew_nbr;

fprintf(['\nSearch time: ', num2str(toc(tSearchStart_nbr), 3), 's\n'])
fprintf('D-opt score at nbr (nominal weights):  %.4f\n', logdet(MoptD_nbr))
fprintf('D-opt score at nbr (optimal weights):  %.4f\n', logdet(MoptD_nbr_true))
fprintf('Efficiency of nominal design at nbr:   %.2f%%\n', ...
        100 * exp((logdet(MoptD_nbr) - logdet(MoptD_nbr_true))/p))


% -------------------------------------------------------------------------
% 6.  Comparison plots
%     Solid lines  = nominal theta   (-o)
%     Dashed lines = neighbour theta (--s)
%     Blue=D, Green=A, Magenta=C, Red=Random
% -------------------------------------------------------------------------

%% Plots
Nexp=100;

    figure(3)   % figure 3 for Nexp=20, figure 7 for Nexp=100
    clf

    % helper: scale by sqrt(Nexp)
    sc = @(M) 1/sqrt(Nexp) * sqrt(diag(inv(M)));
    se = @(M) 1/sqrt(Nexp) * sqrt(eig(inv(M)));

    % --- Subplot 1: uncertainty ellipsoid eigenvalues ---
    subplot(1,3,1)
    semilogy(1:p, se(MoptD),        'b-o',  'LineWidth', 2, 'MarkerSize', 7); hold on
    semilogy(1:p, se(MoptD_nbr),    'b--s', 'LineWidth', 2, 'MarkerSize', 7)
    semilogy(1:p, se(MoptA),        'g-o',  'LineWidth', 2, 'MarkerSize', 7)
    semilogy(1:p, se(MoptA_nbr),    'g--s', 'LineWidth', 2, 'MarkerSize', 7)
    semilogy(1:p, se(MoptC),        'm-o',  'LineWidth', 2, 'MarkerSize', 7)
    semilogy(1:p, se(MoptC_nbr),    'm--s', 'LineWidth', 2, 'MarkerSize', 7)
    semilogy(1:p, se(Mrand_nominal),'r-o',  'LineWidth', 2, 'MarkerSize', 7)
    semilogy(1:p, se(Mrand_nbr),    'r--s', 'LineWidth', 2, 'MarkerSize', 7)
    grid on; xlim([1, p])
    title('Uncertainty ellipsoid')
    ylabel('\surd eigenvalues'); xlabel('Eigenvalue index')
    legend('D nom.','D nbr.','A nom.','A nbr.', ...
           'C nom.','C nbr.','Rand nom.','Rand nbr.','Location','best')
    set(gca, 'FontSize', 10)

    % --- Subplot 2: parameter std dev ---
    subplot(1,3,2)
    plot(1:p, sc(MoptD),        'b-o',  'LineWidth', 2, 'MarkerSize', 7); hold on
    plot(1:p, sc(MoptD_nbr),    'b--s', 'LineWidth', 2, 'MarkerSize', 7)
    plot(1:p, sc(MoptA),        'g-o',  'LineWidth', 2, 'MarkerSize', 7)
    plot(1:p, sc(MoptA_nbr),    'g--s', 'LineWidth', 2, 'MarkerSize', 7)
    plot(1:p, sc(MoptC),        'm-o',  'LineWidth', 2, 'MarkerSize', 7)
    plot(1:p, sc(MoptC_nbr),    'm--s', 'LineWidth', 2, 'MarkerSize', 7)
    plot(1:p, sc(Mrand_nominal),'r-o',  'LineWidth', 2, 'MarkerSize', 7)
    plot(1:p, sc(Mrand_nbr),    'r--s', 'LineWidth', 2, 'MarkerSize', 7)
    grid on; xlim([1, p])
    title('Parameter std dev')
    ylabel('StdDev(\theta_i)'); xlabel('Parameter index')
    legend('D nom.','D nbr.','A nom.','A nbr.', ...
           'C nom.','C nbr.','Rand nom.','Rand nbr.','Location','best')
    set(gca, 'FontSize', 10)
    set(gca, 'XTick', 1:p, 'XTickLabel', param_names, 'TickLabelInterpreter', 'tex')

    % --- Subplot 3: relative std dev (normalised by theta_neighbor) ---
    subplot(1,3,3)
    plot(1:p, sc(MoptD)         ./ thetaData', 'b-o',  'LineWidth', 2, 'MarkerSize', 7); hold on
    plot(1:p, sc(MoptD_nbr)     ./ thetaNeighbor', 'b--s', 'LineWidth', 2, 'MarkerSize', 7)
    plot(1:p, sc(MoptA)         ./ thetaData', 'g-o',  'LineWidth', 2, 'MarkerSize', 7)
    plot(1:p, sc(MoptA_nbr)     ./ thetaNeighbor', 'g--s', 'LineWidth', 2, 'MarkerSize', 7)
    plot(1:p, sc(MoptC)         ./ thetaData', 'm-o',  'LineWidth', 2, 'MarkerSize', 7)
    plot(1:p, sc(MoptC_nbr)     ./ thetaNeighbor', 'm--s', 'LineWidth', 2, 'MarkerSize', 7)
    plot(1:p, sc(Mrand_nominal)  ./ thetaData', 'r-o',  'LineWidth', 2, 'MarkerSize', 7)
    plot(1:p, sc(Mrand_nbr)     ./ thetaNeighbor', 'r--s', 'LineWidth', 2, 'MarkerSize', 7)
    grid on; xlim([1, p])
    title('Relative std dev  (w.r.t. \theta_{nbr})')
    ylabel('StdDev(\theta_i) / \theta_i'); xlabel('Parameter index')
    legend('D nom.','D nbr.','A nom.','A nbr.', ...
           'C nom.','C nbr.','Rand nom.','Rand nbr.','Location','best')
    set(gca, 'FontSize', 10)
    set(gca, 'XTick', 1:p, 'XTickLabel', param_names, 'TickLabelInterpreter', 'tex')

    sgtitle(sprintf('Nominal (solid) vs neighbour (dashed),  \\Delta\\theta = +%.0f%%,  N_{exp} = %d', ...
            perturbationFactor*100, Nexp), 'FontSize', 14, 'FontWeight', 'bold')



    figure(4)   % figure 3 for Nexp=20, figure 7 for Nexp=100
    clf


    % --- Subplot 1: uncertainty ellipsoid eigenvalues ---
    subplot(1,2,1)
    semilogy(1:p, se(MoptD),        'b-o',  'LineWidth', 2, 'MarkerSize', 7); hold on
    semilogy(1:p, se(MoptD_nbr),    'b--s', 'LineWidth', 2, 'MarkerSize', 7)
    semilogy(1:p, se(MoptD_nbr_true),    'b-.d', 'LineWidth', 2, 'MarkerSize', 7)
    % semilogy(1:p, se(MoptA),        'g-o',  'LineWidth', 2, 'MarkerSize', 7)
    % semilogy(1:p, se(MoptA_nbr),    'g--s', 'LineWidth', 2, 'MarkerSize', 7)
    % semilogy(1:p, se(MoptC),        'm-o',  'LineWidth', 2, 'MarkerSize', 7)
    % semilogy(1:p, se(MoptC_nbr),    'm--s', 'LineWidth', 2, 'MarkerSize', 7)
    semilogy(1:p, se(Mrand_nominal),'r-o',  'LineWidth', 2, 'MarkerSize', 7)
    semilogy(1:p, se(Mrand_nbr),    'r--s', 'LineWidth', 2, 'MarkerSize', 7)
    grid on; xlim([1, p])
    title('Uncertainty ellipsoid')
    ylabel('\surd eigenvalues'); xlabel('Eigenvalue index')
    legend('D nom.','D nbr.',...
           'Rand nom.','Rand nbr.','Location','best')
    set(gca, 'FontSize', 10)

    % --- Subplot 2: parameter std dev ---
    subplot(1,2,2)
    plot(1:p, sc(MoptD),        'b-o',  'LineWidth', 2, 'MarkerSize', 7); hold on
    plot(1:p, sc(MoptD_nbr),    'b--s', 'LineWidth', 2, 'MarkerSize', 7)
    plot(1:p, sc(MoptD_nbr_true),    'b-.d', 'LineWidth', 2, 'MarkerSize', 7)
    % plot(1:p, sc(MoptA),        'g-o',  'LineWidth', 2, 'MarkerSize', 7)
    % plot(1:p, sc(MoptA_nbr),    'g--s', 'LineWidth', 2, 'MarkerSize', 7)
    % plot(1:p, sc(MoptC),        'm-o',  'LineWidth', 2, 'MarkerSize', 7)
    % plot(1:p, sc(MoptC_nbr),    'm--s', 'LineWidth', 2, 'MarkerSize', 7)
    plot(1:p, sc(Mrand_nominal),'r-o',  'LineWidth', 2, 'MarkerSize', 7)
    plot(1:p, sc(Mrand_nbr),    'r--s', 'LineWidth', 2, 'MarkerSize', 7)
    grid on; xlim([1, p])
    title('Parameter std dev')
    ylabel('StdDev(\theta_i)'); xlabel('Parameter index')
    legend('D nom.','D nbr.', ...
           'Rand nom.','Rand nbr.','Location','best')
    set(gca, 'FontSize', 10)
    set(gca, 'XTick', 1:p, 'XTickLabel', param_names, 'TickLabelInterpreter', 'tex')
