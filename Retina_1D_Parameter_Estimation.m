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


dn = 50; % Size of one layer

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

delta = retinalWidth/dn; % distance between cells


%%% Determination of size of the noise

% Notice that here, if we keep id as the noise matrix then only the
% relative sizes of the noises matter in the end.
% D is the dynamical noise and Sigmap is the output/measurement noise.

% D = .1 * eye(n); 
% Sigmap = 1 * eye(q); 

% D = 0.5 * eye(n); 
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

% Use a cutoff to accelerate the computations by making the matrix banded:

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

[controls1DWave, inputParameters1DWave] = generate1DPlaneWaveControls(dn, N, Delta, Tf, retinalWidth, cMin, cMax, cN, kMin, kMax, kN);

if addRandomInputs % produce random inputs if asked
    [controlsRandom, inputParametersRandom] = generate1DRandomBinaryControls(dn, N, nbrRandomInputs, seed);
end


% Concatenate inputs and constant input.
if addRandomInputs 
    counter = length(controls1DWave) + length(controlsRandom) + 1;
    inputParameters = [inputParameters1DWave, inputParametersRandom, [NaN;0]];
    controls = cell(counter,1);
    controls(1:end-1) = [controls1DWave; controlsRandom];
else
    counter = length(controls1DWave) + 1;
    inputParameters = [inputParameters1DWave,  [NaN;0]];
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

fprintf(['There are ', num2str(K),' inputs in this experiment. \n'])



%% FISHER MATRICES CONSTRUCTION %%

% Please notice that the Fisher matrices are rescaled 

%%% Collection of Fisher matrices

fprintf(['\n','Fisher matrices computation... \n'])

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

%% OPTIMAL INPUT SELECTION %%

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

wOld = ones(K,1)/K;
wNew = wOld;



% Compute weighted sum
Minit = zeros(p,p);
for k = 1:K
    Minit = Minit + wNew(k) * M{k};
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
    wOld = wNew;

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
    wNew = (1-alpha) * wOld + alpha * wks;
    

    improvement = (logdet(Mnew)-logdet(Mold))/abs(logdet(Mold));

    if mod(i,100) == 0 
        fprintf(['Step of search : ', num2str(i), ', score improvement: ', num2str(improvement), '\n'])
    end
end

tSearchTaken = toc(tSearchStart);

fprintf(['\n','Search time: ', num2str(tSearchTaken,3) ,'s \n'])

fprintf(['Max score of a matrix: ', num2str(maxScore), '\n'])

fprintf(['Final score: ', num2str(logdet(Mnew)), '\n'])

fprintf(['Score of averaged matrix: ', num2str(logdet(Minit)), '\n'])


if addRandomInputs
    wRand = 1/nbrRandomInputs;
    Mrand = zeros(p,p);
    for k = K-1:-1:K-nbrRandomInputs
        Mrand = Mrand + wRand * M{k};
    end
    fprintf(['Score of random inputs: ', num2str(logdet(Mrand)), '\n'])
end


threshold = 0.01;

idx = find(abs(wNew) > threshold);

timeTakenCheck = toc(timeTakenStart);

fprintf(['\n','Total computation time: ', num2str(timeTakenCheck,3) ,'s (',num2str(floor(timeTakenCheck/60)),' min ',num2str(mod(timeTakenCheck,60),2),'s)\n'])



%% REPRESENTATION OF THEORETICAL GAINS %%

Nexp = 100;

figure(1)
clf

% Define parameter names
paramNames = {'\tau_B^{-1}', '\tau_A^{-1}', '\tau_G^{-1}', 'w^-', 'w^+', 'w^B_G','w^A_G'};
paramNamesPlain = {'1/tau_B', '1/tau_A', '1/tau_G', 'w^-', 'w^+', 'w^B_G','w^A_G'};


% --- Subplot 1: Eigenvalues (Comparison) ---
subplot(1,3,1)
semilogy(1:p,1/sqrt(Nexp)*sqrt(eig(inv(Mnew))), 'b-o', 'LineWidth', 2, 'MarkerSize', 8)
hold on
semilogy(1:p, 1/sqrt(Nexp)*sqrt(eig(inv(Mrand))), 'r-s', 'LineWidth', 2, 'MarkerSize', 8)
grid on
title('Predicted Uncertainty ellipsoid')
ylabel('\surd eigenvalues')
xlabel('Eigenvalue Index')
legend('Optimal', 'Random', 'Location', 'best')
set(gca, 'FontSize', 11)

% --- Subplot 2: Parameter Variances (Diagonal of Inverse) ---
subplot(1,3,2)

plot(1:p, 1/sqrt(Nexp)*sqrt(diag(inv(Mnew))), 'b-o', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(1:p, 1/sqrt(Nexp)*sqrt(diag(inv(Mrand))), 'r-s', 'LineWidth', 2, 'MarkerSize', 8)
plot(1:p, 1/sqrt(3)*1/sqrt(Nexp)*sqrt(diag(inv(Mrand))), 'm-.s', 'LineWidth', 2, 'MarkerSize', 8)
grid on
title('Predicted parameter uncertainty')
ylabel('StdDev(\theta_i)')
xlabel('Parameter Index')
legend('Optimal', 'Random', 'Random (rescaled)','Location', 'best')
set(gca, 'FontSize', 11)
set(gca, 'XTick', 1:p, 'XTickLabel', paramNames, 'TickLabelInterpreter', 'tex')

% --- Subplot 3: Relative Parameter Variances ---
subplot(1,3,3)

plot(1:p, 1/sqrt(Nexp)*sqrt(diag(inv(Mnew)))./thetaData', 'b-o', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(1:p, 1/sqrt(Nexp)*sqrt(diag(inv(Mrand)))./thetaData', 'r-s', 'LineWidth', 2, 'MarkerSize', 8)
plot(1:p, 1/sqrt(3)*1/sqrt(Nexp)*sqrt(diag(inv(Mrand)))./thetaData', 'm-.s', 'LineWidth', 2, 'MarkerSize', 8)
grid on
title('Relative predicted parameter uncertainty')
ylabel('StdDev(\theta_i)')
xlabel('Parameter Index')
legend('Optimal', 'Random', 'Random (rescaled)','Location', 'best')
set(gca, 'FontSize', 11)
set(gca, 'XTick', 1:p, 'XTickLabel', paramNames, 'TickLabelInterpreter', 'tex')


sgtitle(sprintf('Predicted uncertainty for D-Optimal vs Random Inputs  (K=%d)', Nexp), ...
        'FontSize', 16, 'FontWeight', 'bold')



%% PLOTTING SOME TRAJECTORIES %%

figure(2)
clf

Nexp=100;

[wOptimal, sortOrder] = sort(wNew(idx), 'descend');
idxOptimal = idx(sortOrder);
nnOptimal = efficientRounding(wOptimal, Nexp);
nnUniform = ones(1,4);

SigmaBold = buildSigmaBold(C, C/2, F, Sigma, Sigmap, N);

rng(312)
[Ybar, Eta, SigmaBoldHat] = dataStochastic1D(controls, idxOptimal, nnOptimal, q, C, F, G, SigmaBold, N);
[YbarUnif, EtaUnif, SigmaBoldHatUnif] = dataStochastic1D(controls, idxOptimal, nnUniform, q, C, F, G, SigmaBold, N);

IntT = (1:N)*Delta/1000;


yLimInf = 10+ceil(max(abs(YbarUnif(1:q:end,:)),[],"all")/10)*10;

inputPosition = zeros(4,N);
for i = 1:min(4,length(idxOptimal))
    inputPosition(i,:) = controls{idxOptimal(i)}(positions(1),:);
    inputPosition(i,:) = inputPosition(i,:)/max(inputPosition(i,:));
end

for i = 1:min(4,length(nnOptimal))
    subplot(2,2,i)
    
    plot(IntT,20*inputPosition(i,:),'k-', 'LineWidth', 1)
    hold on
    plot(IntT,Eta(1:q:end,i),'k:', 'LineWidth', 1)
    plot(IntT,YbarUnif(1:q:end,i),'k-.', 'LineWidth', 1)
    plot(IntT,Ybar(1:q:end,i),'r-', 'LineWidth', 1)
    ylim([-yLimInf,30])
    xlim([0,IntT(end)])
    title(sprintf('Weight = %0.3f, %2d repetitions', [wOptimal(i),nnOptimal(i)]))
    ylabel('Output')
    xlabel('Time (sec)')
    legend('Input (rescaled)', 'Deterministic response', 'Output', 'Output avg over repet.', 'Location', 'best')
    set(gca, 'FontSize', 11)
end
sgtitle('Few input and output evolutions (D-optimal case)', 'FontSize', 16, 'FontWeight', 'bold')


figure(3)
clf

idxRdm = K-1:-1:K-1-Nexp+1;
nnRdm = ones(1, length(idxRdm));

[Ybar, Eta, SigmaBoldHat] = dataStochastic1D(controls, idxRdm, nnRdm, q, C, F, G, SigmaBold, N);

for i = 1:4
    inputPosition(i,:) = controls{idxRdm(i)}(positions(1),:);
    inputPosition(i,:) = inputPosition(i,:)/max(inputPosition(i,:));
end

%yLimInf = 10+ceil(max(abs(Ybar(1:q:end,:)),[],"all")/10)*10;

for i = 1:4
    subplot(2,2,i)
    stairs(IntT,20*inputPosition(i,:),'k-', 'LineWidth', 1)
    hold on
    plot(IntT(1:end-1),Eta(q+1:q:end,i),'k:', 'LineWidth', 1)
    plot(IntT,Ybar(1:q:end,i),'k-.', 'LineWidth', 1)
    ylim([-yLimInf,30])
    xlim([0,IntT(end)])
    title(' ')
    ylabel('Output')
    xlabel('Time (sec)')
    legend('Input, rescaled', 'Deterministic response', 'Output', 'Location', 'best')
    set(gca, 'FontSize', 11)
end
sgtitle('Few input and output evolutions (Random case)', 'FontSize', 16, 'FontWeight', 'bold')


%% ESTIMATION PREAMBLE %%

fprintf('\n========================================\n')
fprintf('Comparative study: D-optimal vs random inputs.\n')
fprintf('========================================\n\n')



Nexp = 100;
Nrepeat = 150;  % Increase for better statistics

thetaInit = ones(p,1)*0.01;


if cutoff
    SigmaBold = buildSigmaBoldCutoff(C, C/2, F, Sigma, Sigmap, N, depth);
else
    SigmaBold = buildSigmaBold(C, C/2, F, Sigma, Sigmap, N);
end


% D-optimal inputs (already defined from optimization)

wApproximate = wNew(idx);
idxOptimal = idx;
nnOptimal = efficientRounding(wApproximate, Nexp);

% Random inputs

idxRdm = [counter-Nexp-1:1:counter-1];
nnRdm = ones(1, Nexp);

% Store results
Results = struct();

% Optimization settings for fmincon

displayOption = 'iter-detailed';
displayOption = 'final-detailed';
algorithmOption = 'interior-point';

options = optimoptions(@fmincon, ...
    'Algorithm', algorithmOption, ...
    'MaxIterations', 1000, ...          
    'MaxFunctionEvaluations', 5000, ...
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-10, ...
    'Display', displayOption, ...
    'ConstraintTolerance', 1e-6);


%% D-OPTIMAL INPUTS %%
fprintf(['\n','D-optimal inputs... \n'])


thetaOptimalAll = NaN(Nrepeat, p);
exitFlagsOptimal = NaN(Nrepeat, 1);
fValsOptimal = NaN(Nrepeat, 1);


timebar = waitbar(0, 'D-optimal: repeated simulations...');
tLoopStart = tic;

for i = 1:Nrepeat
    
    % Generate stochastic data for this iteration
    [Ybar, Eta, SigmaBoldHat] = dataStochastic1D(controls, idxOptimal, nnOptimal, q, C, F, G, SigmaBold, N);
    
    % Objective function    
    JML = @(theta) JEfficientStats(theta, Ybar, SigmaBoldHat, controls, idxOptimal, nnOptimal, Delta, Acell, B, C, N, D, Sigmap);
    
    % Box constraints: 0.5x and 1.5x of the true thetaData
    lb = thetaData' * 0.5; 
    ub = thetaData' * 1.5;

    [thetaOpt, fVal, exitFlag] = fmincon(JML, thetaInit, [], [], [], [], lb, ub, [], options);
    
    thetaOptimalAll(i, :) = thetaOpt;
    exitFlagsOptimal(i) = exitFlag;
    fValsOptimal(i) = fVal;

    timeInLoop = toc(tLoopStart);
    timeRemainingMin = floor((Nrepeat-i)/i*timeInLoop/60);
    timeRemainingSec = floor((Nrepeat-i)/i*timeInLoop-60*timeRemainingMin);
    
    waitbar(i/Nrepeat, timebar, ['D-optimal (',...
    num2str(floor(i/Nrepeat*100)),' %, est. remaining time: ',...
    num2str(timeRemainingMin),'min ',...
    num2str(timeRemainingSec),'s)']);
end


timeOptimal = toc(tLoopStart);
close(timebar);

% Compute statistics
[meanThetaOpt, meanAbsBiasOpt, covThetaOpt, sqrtEigCovOpt] = ...
    statTheta(thetaOptimalAll, thetaData);

Results.optimal.meanTheta = meanThetaOpt;
Results.optimal.meanAbsBias = meanAbsBiasOpt;
Results.optimal.covTheta = covThetaOpt;
Results.optimal.sqrtEigCov = sqrtEigCovOpt;
Results.optimal.exitFlags = exitFlagsOptimal;
Results.optimal.computationTime = timeOptimal;

fprintf(['D-optimal completed in ', num2str(timeOptimal,3), 's. \n\n'])


%% RANDOM INPUTS %%

fprintf(['\n','Random inputs... \n'])


thetaRandomAll = NaN(Nrepeat, p);
exitFlagsRandom = NaN(Nrepeat, 1);
fValsRandom = NaN(Nrepeat, 1);


timebar = waitbar(0, 'Random: repeated simulations...');
tLoopStart = tic;

for i = 1:Nrepeat
    
    % Generate stochastic data for this iteration
    [Ybar, Eta, SigmaBoldHat] = dataStochastic1D(controls, idxRdm, nnRdm, q, C, F, G, SigmaBold, N);
    
    % Objective function    
    JML = @(theta) JEfficientStats(theta, Ybar, SigmaBoldHat, controls, idxRdm, nnRdm, Delta, Acell, B, C, N, D, Sigmap);
    
    % Box constraints: 0.5x and 1.5x of the true thetaData
    lb = thetaData' * 0.5; 
    ub = thetaData' * 1.5;

    [thetaOpt, fVal, exitFlag] = fmincon(JML, thetaInit, [], [], [], [], lb, ub, [], options);
    
    thetaRandomAll(i, :) = thetaOpt;
    exitFlagsRandom(i) = exitFlag;
    fValsRandom(i) = fVal;

    timeInLoop = toc(tLoopStart);
    timeRemainingMin = floor((Nrepeat-i)/i*timeInLoop/60);
    timeRemainingSec = floor((Nrepeat-i)/i*timeInLoop-60*timeRemainingMin);
    
    waitbar(i/Nrepeat, timebar, ['Random (',...
    num2str(floor(i/Nrepeat*100)),' %, est. remaining time: ',...
    num2str(timeRemainingMin),'min ',...
    num2str(timeRemainingSec),'s)']);
    
end


timeRandom = toc(tLoopStart);
close(timebar);

% Compute statistics
[meanThetaRdm, meanAbsBiasRdm, covThetaRdm, sqrtEigCovRdm] = ...
    statTheta(thetaRandomAll, thetaData);

Results.random.meanTheta = meanThetaRdm;
Results.random.meanAbsBias = meanAbsBiasRdm;
Results.random.covTheta = covThetaRdm;
Results.random.sqrtEigCov = sqrtEigCovRdm;
Results.random.exitFlags = exitFlagsRandom;
Results.random.computationTime = timeRandom;

fprintf(['Random inputs completed in ', num2str(timeRandom,3), 's. \n\n'])

%% DISPLAY RESULTS %%

fprintf(['\n','========================================\n'])
fprintf(['Results summary.\n'])
fprintf(['========================================\n\n'])


% Mean absolute bias
fprintf(['Mean absolute bias (averaged across parameters): \n'])
fprintf(['  D-optimal: ', num2str(mean(meanAbsBiasOpt),3), '\n'])
fprintf(['  Random:    ', num2str(mean(meanAbsBiasRdm),3), '\n'])
fprintf(['  Ratio (Random/Optimal): ', num2str(mean(meanAbsBiasRdm)/mean(meanAbsBiasOpt),3), '\n\n'])

% Uncertainty (geometric mean of eigenvalues)
geomeanOpt = exp(mean(log(sqrtEigCovOpt)));
geomeanRdm = exp(mean(log(sqrtEigCovRdm)));

fprintf(['Uncertainty (geometric mean of sqrt eigenvalues): \n'])
fprintf(['  D-optimal: ', num2str(geomeanOpt,3), '\n'])
fprintf(['  Random:    ', num2str(geomeanRdm,3), '\n'])
fprintf(['  Ratio (Random/Optimal): ', num2str(geomeanRdm/geomeanOpt,3), '\n\n'])

% Determinant of covariance (volume of uncertainty ellipsoid)
detOpt = det(covThetaOpt);
detRdm = det(covThetaRdm);

fprintf(['Covariance determinant (uncertainty volume): \n'])
fprintf(['  D-optimal: ', num2str(detOpt,3), '\n'])
fprintf(['  Random:    ', num2str(detRdm,3), '\n'])
fprintf(['  Ratio (Random/Optimal): ', num2str(detRdm/detOpt,3), '\n\n'])

% Parameter-by-parameter comparison
fprintf(['Parameter-by-parameter standard deviations: \n'])
fprintf(['  Param.     D-optimal    Random       Ratio (R/O)\n'])
fprintf(['  ---------------------------------------------------\n'])
for i = 1:p
    stdOpt = sqrt(covThetaOpt(i,i));
    stdRdm = sqrt(covThetaRdm(i,i));
    fprintf(['  θ', num2str(i), '         ', num2str(stdOpt,3), '    ', num2str(stdRdm,3), '    ', num2str(stdRdm/stdOpt,3), '\n'])
end
fprintf(['\n'])

% Eigenvalue comparison
fprintf(['Eigenvalues comparison (std devs along principal axes): \n'])
fprintf(['  Axis       D-optimal    Random       Ratio (R/O)\n'])
fprintf(['  ---------------------------------------------------\n'])
for i = 1:p
    fprintf(['  ', num2str(i), '           ', num2str(sqrtEigCovOpt(i),3), '    ', num2str(sqrtEigCovRdm(i),3), '    ', num2str(sqrtEigCovRdm(i)/sqrtEigCovOpt(i),3), '\n'])
end
fprintf(['\n'])



%%  VISUALIZATIONS %%
figure(4)
clf

% 1. Confidence Intervals Comparison (Normalized)
subplot(2,3,1)

ciWidth = 1.96; % for 95% CI

stdOpt = sqrt(diag(covThetaOpt));
stdRdm = sqrt(diag(covThetaRdm));

% Normalize by true values
meanNormOpt = meanThetaOpt ./ thetaData;
meanNormRdm = meanThetaRdm ./ thetaData;

stdNormOpt = stdOpt ./ thetaData';
stdNormRdm = stdRdm ./ thetaData';

ciLowerOpt = meanNormOpt' - ciWidth * stdNormOpt;
ciUpperOpt = meanNormOpt' + ciWidth * stdNormOpt;

ciLowerRdm = meanNormRdm' - ciWidth * stdNormRdm;
ciUpperRdm = meanNormRdm' + ciWidth * stdNormRdm;

hold on
xPositions = 1:p;
offset = 0.2;

barWidth = 0.1;

for i = 1:p
    % D-optimal (left, blue)
    plot([i-offset, i-offset], [ciLowerOpt(i), ciUpperOpt(i)], ...
         'b-', 'LineWidth', 1.5)
    plot(i-offset, meanNormOpt(i), 'b+', 'MarkerSize', 10, 'LineWidth', 1.5)
    plot([i-offset-barWidth, i-offset+barWidth], [ciLowerOpt(i), ciLowerOpt(i)], ...
         'b-', 'LineWidth', 1.5)
    plot([i-offset-barWidth, i-offset+barWidth], [ciUpperOpt(i), ciUpperOpt(i)], ...
         'b-', 'LineWidth', 1.5)

    % Random (right, red)
    plot([i+offset, i+offset], [ciLowerRdm(i), ciUpperRdm(i)], ...
         'r-', 'LineWidth', 1.5)
    plot(i+offset, meanNormRdm(i), 'r+', 'MarkerSize', 10, 'LineWidth', 1.5)
    plot([i+offset-barWidth, i+offset+barWidth], [ciLowerRdm(i), ciLowerRdm(i)], ...
         'r-', 'LineWidth', 1.5)
    plot([i+offset-barWidth, i+offset+barWidth], [ciUpperRdm(i), ciUpperRdm(i)], ...
         'r-', 'LineWidth', 1.5)
end

set(gca, 'XTick', 1:p, 'XTickLabel', paramNames, 'TickLabelInterpreter', 'tex')
ylabel('Normalized Parameter Value')
title('95% Confidence Intervals (Normalized)')
yline(1, 'k--', 'LineWidth', 1, 'Alpha', 0.5)
legend('', 'D-optimal mean', '', '', '', 'Random mean', '', '', ...
       'Location', 'best')
grid on
set(gca, 'FontSize', 11)
hold off


% 2. Bias comparison
subplot(2,3,2)
bar([meanAbsBiasOpt; meanAbsBiasRdm]')
set(gca, 'XTickLabel', paramNames, 'TickLabelInterpreter', 'tex')
ylabel('Mean Absolute Bias')
title('Parameter Estimation Bias')
legend('D-optimal', 'Random', 'Location', 'best')
grid on
set(gca, 'FontSize', 11)

% 3. Standard deviations comparison
subplot(2,3,3)
stdOpt = sqrt(diag(covThetaOpt));
stdRdm = sqrt(diag(covThetaRdm));
bar([stdOpt, stdRdm])
set(gca, 'XTickLabel', paramNames, 'TickLabelInterpreter', 'tex')
ylabel('Standard Deviation')
title('Parameter Uncertainty')
legend('D-optimal', 'Random', 'Location', 'best')
grid on
set(gca, 'FontSize', 11)

% 4. Eigenvalue comparison
subplot(2,3,4)
semilogy(1:p, sqrtEigCovOpt, 'b-o', 'LineWidth', 2, 'MarkerSize', 8)
hold on
semilogy(1:p, sqrtEigCovRdm, 'r-s', 'LineWidth', 2, 'MarkerSize', 8)
xlabel('Eigenvalue Index')
ylabel('\surd(Eigenvalue) [log scale]', 'Interpreter', 'tex')
title('Uncertainty Ellipsoid Principal Axes')
legend('D-optimal', 'Random', 'Location', 'best')
grid on
set(gca, 'FontSize', 11)

% 5. Parameter estimates scatter (first 2 parameters)
subplot(2,3,5)
plot(thetaOptimalAll(:,1), thetaOptimalAll(:,2), 'b.', 'MarkerSize', 10)
hold on
plot(thetaRandomAll(:,1), thetaRandomAll(:,2), 'r.', 'MarkerSize', 10)
plot(thetaData(1), thetaData(2), 'k*', 'MarkerSize', 15, 'LineWidth', 2)
xlabel(paramNames{1}, 'Interpreter', 'tex', 'FontSize', 12)
ylabel(paramNames{2}, 'Interpreter', 'tex', 'FontSize', 12)
title(['Parameter Estimates (', paramNames{1}, ' vs ', paramNames{2}, ')'], ...
      'Interpreter', 'tex')
legend('D-optimal', 'Random', 'True value', 'Location', 'best')
grid on
set(gca, 'FontSize', 11)

% 6. Ratio of uncertainties
subplot(2,3,6)
ratios = sqrtEigCovRdm ./ sqrtEigCovOpt;
bar(ratios)
xlabel('Eigenvalue Index')
ylabel('Ratio (Random/D-optimal)')
title('Relative Efficiency')
yline(1, 'k--', 'LineWidth', 1.5)
grid on
set(gca, 'FontSize', 11)

sgtitle(sprintf('D-Optimal vs Random Inputs Comparison (N=%d repeats)', Nrepeat), ...
        'FontSize', 16, 'FontWeight', 'bold')


%% PARAMETER-SPECIFIC DISTRIBUTIONS %%

figure(5)
clf

for i = 1:p
    subplot(2, 4, i)
    
    histogram(thetaOptimalAll(:,i), 15, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'Normalization', 'pdf')
    hold on
    histogram(thetaRandomAll(:,i), 15, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'Normalization', 'pdf')
    
    xline(thetaData(i), 'k--', 'LineWidth', 2, 'Label', 'True', ...
          'LabelHorizontalAlignment', 'center', 'FontSize', 9)
    
    xlabel(['Parameter ', paramNames{i}], 'Interpreter', 'tex', 'FontSize', 11)
    ylabel('Probability Density')
    title(paramNames{i}, 'Interpreter', 'tex', 'FontSize', 13, 'FontWeight', 'bold')
    
    if i == 1
        legend('D-optimal', 'Random', 'Location', 'best')
    end
    
    grid on
    set(gca, 'FontSize', 10)
end

sgtitle('Distribution of Parameter Estimates', 'FontSize', 16, 'FontWeight', 'bold')


%% PREDICTED VS ESTIMATED UNCERTAINTY %%

figure(6)
clf

% --- Subplot 1: Eigenvalues (Comparison) ---
subplot(1,2,1)
semilogy(1:p, 1/sqrt(Nexp)*sqrt(eig(inv(Mnew))), 'b--o', 'LineWidth', 2, 'MarkerSize', 8)
hold on
semilogy(1:p, 1/sqrt(Nexp)*sqrt(eig(inv(Mrand))), 'r--s', 'LineWidth', 2, 'MarkerSize', 8)
semilogy(1:p, sqrtEigCovOpt, 'b-o', 'LineWidth', 2, 'MarkerSize', 8)
semilogy(1:p, sqrtEigCovRdm, 'r-s', 'LineWidth', 2, 'MarkerSize', 8)


grid on
title('Predicted vs estimated uncertainty ellipsoid')
ylabel('\surd eigenvalues')
xlabel('Eigenvalue Index')
legend('D-Optimal (predicted)', 'Random (predicted)', 'D-Optimal (estimated)', 'Random (estimated)', 'Location', 'northwest')
set(gca, 'FontSize', 11)

% --- Subplot 2: Parameter Variances (Diagonal of Inverse) ---
subplot(1,2,2)

plot(1:p, 1/sqrt(Nexp)*sqrt(diag(inv(Mnew))), 'b--o', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(1:p, 1/sqrt(Nexp)*sqrt(diag(inv(Mrand))), 'r--s', 'LineWidth', 2, 'MarkerSize', 8)
plot(1:p, sqrt(diag(covThetaOpt)), 'b-o', 'LineWidth', 2, 'MarkerSize', 8)
plot(1:p, sqrt(diag(covThetaRdm)), 'r-s', 'LineWidth', 2, 'MarkerSize', 8)
plot(1:p, 1/sqrt(3)*1/sqrt(Nexp)*sqrt(diag(inv(Mrand))), 'r:s', 'LineWidth', 2, 'MarkerSize', 8)
grid on
title('Predicted vs estimated parameter uncertainty')
ylabel('StdDev(\theta_i)')
xlabel('Parameter Index')
legend('D-Optimal (predicted)', 'Random (predicted)', 'D-Optimal (estimated)', 'Random (estimated)','Random rescaled (predicted)', 'Location', 'northeast')
set(gca, 'FontSize', 11)
set(gca, 'XTick', 1:p, 'XTickLabel', paramNames, 'TickLabelInterpreter', 'tex')




%% PAIRWISE PARAMETER SCATTER WITH CONFIDENCE ELLIPSES %%

figure(7)
clf
%set(figure(7), 'Position', [100, 100, 1400, 900])

nPairs = p*(p-1)/2;
nCols  = 6;
nRows  = ceil(nPairs/nCols);

chiSq   = -2 * log(0.05);
padding = 1.5;

% Per-parameter half-width in original units
halfWidth = padding * sqrt(chiSq) * sqrt(diag(covThetaRdm))';  % (1 x p)

iPair = 0;

for i = 1:p
    for j = i+1:p

        iPair = iPair + 1;
        subplot(nRows, nCols, iPair)
        hold on

        % Point clouds in original units
        plot(thetaOptimalAll(:,i), thetaOptimalAll(:,j), 'b.', 'MarkerSize', 4)
        plot(thetaRandomAll(:,i),  thetaRandomAll(:,j),  'r.', 'MarkerSize', 4)

        

        % Confidence ellipses in original units
        plotConfidenceEllipse(covThetaOpt([i,j],[i,j]), ...
            [meanThetaOpt(i), meanThetaOpt(j)], 0.95, 'b')
        plotConfidenceEllipse(covThetaRdm([i,j],[i,j]), ...
            [meanThetaRdm(i), meanThetaRdm(j)], 0.95, 'r')

        % True value
        plot(thetaData(i), thetaData(j), 'w+', 'MarkerSize', 10, 'LineWidth', 2)
        plot(thetaData(i), thetaData(j), 'k+', 'MarkerSize', 9, 'LineWidth', 1.)

        % Use the LARGER of the two half-widths for both axes,
        % so that unit lengths are equal on x and y within each panel.
        hw = max(halfWidth(i), halfWidth(j));
        xlim(thetaData(i) + hw * [-1, 1])
        ylim(thetaData(j) + hw * [-1, 1])

        axis equal
        box on

        xlabel(paramNames{i}, 'Interpreter', 'tex', 'FontSize', 9)
        ylabel(paramNames{j}, 'Interpreter', 'tex', 'FontSize', 9)
        grid on
        set(gca, 'FontSize', 8)

        % if iPair == 1
        %     legend('D-optimal', 'Random', 'True value', ...
        %         'Location', 'best', 'FontSize', 7)
        % end

    end
end

sgtitle('Pairwise parameter estimates and 95% confidence ellipses', ...
    'FontSize', 14, 'FontWeight', 'bold')
%% SAVE RESULTS %%

% Create timestamped subfolder
timeStamp = datestr(now, 'yyyy-mm-dd-HH-MM');
saveDir   = fullfile('EstimationResults', timeStamp);
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

fprintf(['\n', 'Saving results to: ', saveDir, '... \n'])

% --- File 1: Estimation results ---
% Everything needed to analyse and re-plot the comparison.
estimationResults = struct();

estimationResults.thetaOptimalAll   = thetaOptimalAll;
estimationResults.thetaRandomAll    = thetaRandomAll;
estimationResults.exitFlagsOptimal  = exitFlagsOptimal;
estimationResults.exitFlagsRandom   = exitFlagsRandom;
estimationResults.fValsOptimal      = fValsOptimal;
estimationResults.fValsRandom       = fValsRandom;
estimationResults.Results           = Results;

estimationResults.Nrepeat           = Nrepeat;
estimationResults.Nexp              = Nexp;
estimationResults.idxOptimal        = idxOptimal;
estimationResults.nnOptimal         = nnOptimal;
estimationResults.idxRdm            = idxRdm;
estimationResults.nnRdm             = nnRdm;
estimationResults.wNew              = wNew;
estimationResults.idx               = idx;

estimationResults.logDetOptimal     = logdet(Mnew);
estimationResults.logDetRandom      = logdet(Mrand);
estimationResults.logDetInit        = logdet(Minit);
estimationResults.Mnew              = Mnew;
estimationResults.Mrand             = Mrand;
estimationResults.Minit             = Minit;

save(fullfile(saveDir, 'estimationResults.mat'), '-struct', 'estimationResults');
fprintf(['  estimationResults.mat saved. \n'])

% --- File 2: Model and experiment configuration ---
% Everything needed to reproduce the experiment exactly from scratch.
modelConfig = struct();

% True parameter vector and names
modelConfig.thetaData           = thetaData;
modelConfig.thetaNominal        = thetaNominal;
modelConfig.paramNames          = paramNames;
modelConfig.paramNamesPlain     = paramNamesPlain;

% System matrices
modelConfig.A                   = full(A);
modelConfig.B                   = full(B);
modelConfig.C                   = full(C);
modelConfig.C0                  = full(C0);
modelConfig.D                   = D;
modelConfig.Sigmap              = Sigmap;
modelConfig.F                   = full(F);
modelConfig.G                   = full(G);
modelConfig.Sigma               = Sigma;

% Topology matrices
modelConfig.Gamma               = full(Gamma);
modelConfig.GammaBG             = full(GammaBG);
modelConfig.GammaAG             = full(GammaAG);

% Acell (needed to rerun JEfficientStats)
modelConfig.Acell               = cellfun(@full, Acell, 'UniformOutput', false);

% Physical parameters
modelConfig.tauB                = tauB;
modelConfig.tauA                = tauA;
modelConfig.tauG                = tauG;
modelConfig.wm                  = wm;
modelConfig.wp                  = wp;
modelConfig.wBG                 = wBG;
modelConfig.wAG                 = wAG;
modelConfig.retinalWidth        = retinalWidth;
modelConfig.sigmaBG             = sigmaBG;
modelConfig.sigmaAG             = sigmaAG;

% Numerical parameters
modelConfig.dn                  = dn;
modelConfig.n                   = n;
modelConfig.m                   = m;
modelConfig.p                   = p;
modelConfig.q                   = q;
modelConfig.N                   = N;
modelConfig.Delta               = Delta;
modelConfig.DeltaMin            = DeltaMin;
modelConfig.Tf                  = Tf;
modelConfig.positions           = positions;
modelConfig.delta               = delta;

% RNG and reproducibility
modelConfig.seed                = seed;
modelConfig.cutoff              = cutoff;
if cutoff
    modelConfig.depth           = depth;
end
modelConfig.addRandomInputs     = addRandomInputs;
modelConfig.nbrRandomInputs     = nbrRandomInputs;

% Optimisation settings
modelConfig.algorithmOption     = algorithmOption;
modelConfig.displayOption       = 'none';   % don't save 'iter-detailed'
modelConfig.thetaInit           = thetaInit;

% Timestamp and environment
modelConfig.timeStamp           = timeStamp;
modelConfig.matlabVersion       = version;

save(fullfile(saveDir, 'modelConfig.mat'), '-struct', 'modelConfig');
fprintf(['  modelConfig.mat saved. \n'])

% --- File 3: Human-readable summary ---
summaryFile = fopen(fullfile(saveDir, 'summary.txt'), 'w');

fprintf(summaryFile, '========================================\n');
fprintf(summaryFile, 'Experiment summary\n');
fprintf(summaryFile, 'Timestamp: %s\n', timeStamp);
fprintf(summaryFile, 'MATLAB version: %s\n', version);
fprintf(summaryFile, '========================================\n\n');

fprintf(summaryFile, 'MODEL\n');
fprintf(summaryFile, '  State dimension n     : %d (= 3 x dn, dn = %d)\n', n, dn);
fprintf(summaryFile, '  Output dimension q    : %d\n', q);
fprintf(summaryFile, '  Input dimension m     : %d\n', m);
fprintf(summaryFile, '  Parameters p          : %d\n', p);
fprintf(summaryFile, '  Time steps N          : %d\n', N);
fprintf(summaryFile, '  Sampling period Delta : %.4f ms\n', Delta);
fprintf(summaryFile, '  Total duration Tf     : %.2f ms\n', Tf);
fprintf(summaryFile, '  Retinal width         : %.2f mm\n', retinalWidth);
fprintf(summaryFile, '  Output positions      : %s\n', mat2str(positions));
fprintf(summaryFile, '\n');

fprintf(summaryFile, 'TRUE PARAMETERS\n');
for i = 1:p
    fprintf(summaryFile, '  %-12s = %.8f\n', paramNamesPlain{i}, thetaData(i));
end
fprintf(summaryFile, '\n');

fprintf(summaryFile, 'NOISE\n');
fprintf(summaryFile, '  D      : %.4f * I_%d\n', D(1,1), n);
fprintf(summaryFile, '  Sigmap : %.4f * I_%d\n', Sigmap(1,1), q);
fprintf(summaryFile, '\n');

fprintf(summaryFile, 'EXPERIMENT DESIGN\n');
fprintf(summaryFile, '  Nexp (inputs per experiment) : %d\n', Nexp);
fprintf(summaryFile, '  Nrepeat                      : %d\n', Nrepeat);
fprintf(summaryFile, '  RNG seed                     : %d\n', seed);
fprintf(summaryFile, '  Cutoff                       : %d\n', cutoff);
fprintf(summaryFile, '  Random inputs added          : %d\n', addRandomInputs);
if addRandomInputs
    fprintf(summaryFile, '  Number of random inputs      : %d\n', nbrRandomInputs);
end
fprintf(summaryFile, '\n');

fprintf(summaryFile, 'OPTIMAL DESIGN\n');
fprintf(summaryFile, '  log det(M_optimal) : %.6f\n', logdet(Mnew));
fprintf(summaryFile, '  log det(M_random)  : %.6f\n', logdet(Mrand));
fprintf(summaryFile, '  log det(M_init)    : %.6f\n', logdet(Minit));
fprintf(summaryFile, '  Active inputs (weight > %.2f):\n', threshold);
for i = 1:length(idxOptimal)
    fprintf(summaryFile, '    Input %d: weight = %.4f, repetitions = %d\n', ...
        idxOptimal(i), wNew(idxOptimal(i)), nnOptimal(i));
end
fprintf(summaryFile, '\n');

fprintf(summaryFile, 'RESULTS SUMMARY\n');
successOpt = sum(exitFlagsOptimal == 1);
successRdm = sum(exitFlagsRandom == 1);
fprintf(summaryFile, '  Optimisation success rates:\n');
fprintf(summaryFile, '    D-optimal : %d/%d (%.1f%%)\n', successOpt, Nrepeat, 100*successOpt/Nrepeat);
fprintf(summaryFile, '    Random    : %d/%d (%.1f%%)\n', successRdm, Nrepeat, 100*successRdm/Nrepeat);
fprintf(summaryFile, '\n');
fprintf(summaryFile, '  Parameter-by-parameter standard deviations:\n');
fprintf(summaryFile, '  %-12s  %-12s  %-12s  %-10s\n', 'Parameter', 'D-optimal', 'Random', 'Ratio R/O');
fprintf(summaryFile, '  %s\n', repmat('-', 1, 52));
for i = 1:p
    stdOpt = sqrt(Results.optimal.covTheta(i,i));
    stdRdm = sqrt(Results.random.covTheta(i,i));
    fprintf(summaryFile, '  %-12s  %-12.4e  %-12.4e  %-10.3f\n', ...
        paramNamesPlain{i}, stdOpt, stdRdm, stdRdm/stdOpt);
end
fprintf(summaryFile, '\n');
fprintf(summaryFile, '  Geometric mean uncertainty:\n');
fprintf(summaryFile, '    D-optimal : %.4e\n', exp(mean(log(Results.optimal.sqrtEigCov))));
fprintf(summaryFile, '    Random    : %.4e\n', exp(mean(log(Results.random.sqrtEigCov))));
fprintf(summaryFile, '\n');
fprintf(summaryFile, '  Computation times:\n');
fprintf(summaryFile, '    D-optimal : %.1f s\n', Results.optimal.computationTime);
fprintf(summaryFile, '    Random    : %.1f s\n', Results.random.computationTime);

fclose(summaryFile);
fprintf(['  summary.txt saved. \n'])

% --- Figures ---
% Figures as they appear in the script:
%   1 - Predicted uncertainty (D-optimal vs random, 3 subplots)
%   2 - Optimal input trajectories
%   3 - Random input trajectories
%   4 - Comparison results (6 subplots: CIs, bias, std devs, eigenvalues, scatter, efficiency)
%   5 - Parameter estimate distributions (histograms)
%   6 - Predicted vs estimated uncertainty
%   7 - Pairwise scatter plots 

figureNames = { ...
    'fig1_predictedUncertainty', ...
    'fig2_optimalTrajectories', ...
    'fig3_randomTrajectories', ...
    'fig4_comparisonResults', ...
    'fig5_parameterDistributions', ...
    'fig6_predictedVsEstimated',...
    'fig7_scatterPlots'};

for iFig = 1:7
    fig = figure(iFig);
    if ~isempty(fig) && isvalid(fig)
        savefig(fig, fullfile(saveDir, [figureNames{iFig}, '.fig']));
        exportgraphics(fig, fullfile(saveDir, [figureNames{iFig}, '.pdf']), ...
            'ContentType', 'vector');
    end
end
fprintf(['  Figures 1-7 saved. \n'])

fprintf(['\n', 'All results saved to: ', saveDir, '. \n'])
fprintf(['\n', '========================================\n'])
fprintf(['Comparison complete. \n'])
fprintf(['========================================\n'])