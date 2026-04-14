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

delta = retinalWidth/dn; % distance between cells


%%% Determination of size of the noise

% Notice that here, if we keep id as the noise matrix then only the
% relative sizes of the noises matter in the end.
% D is the dynamical noise and Sigmap is the output/measurement noise.

% D = .1 * eye(n); 
% Sigmap = 1 * eye(q); 

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

cN = 200;
kN = 200;

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



%% representation of theoretical gains

Nexp = 20;

figure(1)
clf


% Define parameter names
param_names = {'\tau_B^{-1}', '\tau_A^{-1}', '\tau_G^{-1}', 'w^-', 'w^+', 'w^B_G','w^A_G'};
param_names_plain = {'1/tau_B', '1/tau_A', '1/tau_G', 'w^-', 'w^+', 'w^B_G','w^A_G'};


% --- Subplot 1: Eigenvalues (Comparison) ---
subplot(1,3,1)
% Use semilogy to match the style of your other eigenvalue plot
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

% Plotting diag(inv(M)) gives the variance of each parameter
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
set(gca, 'XTick', 1:p, 'XTickLabel', param_names, 'TickLabelInterpreter', 'tex')

% --- Subplot 2: Parameter Variances (Diagonal of Inverse) ---
subplot(1,3,3)

% Plotting diag(inv(M)) gives the variance of each parameter
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
set(gca, 'XTick', 1:p, 'XTickLabel', param_names, 'TickLabelInterpreter', 'tex')


sgtitle(sprintf('Predicted uncertainty for D-Optimal vs Random Inputs  (K=%d)', Nexp), ...
        'FontSize', 16, 'FontWeight', 'bold')




%% ANIMATE THE INPUTS %%

figure(1)
clf
pause(.5)

% Define spatial domain

kStar = 2*pi/retinalWidth;

dx = retinalWidth/(dn-1);
x = 0:dx:(dn-1)*dx;

lidx = length(idx);
info = cell(lidx,3);
for k=1:lidx
    info{k,1} = num2str(wnew(idx(k)));
    info{k,2} = num2str(inputParameters(1,idx(k))/cStar);
    info{k,3} = num2str(inputParameters(2,idx(k))/kStar);
end

control_norm = zeros(lidx,1);
for k = 1:lidx
    fprintf("weight = " + info{k,1} + ...
            ", c = " + info{k,2} + ...
            " x c^*, k = " + info{k,3} + " x k^*.\n")
    control_norm(k) = max(controls{idx(k)},[],"all");
end




for l=1:N
    for k=1:lidx
        subplot(length(idx),1,k)
        control_plot = controls{idx(k)}/control_norm(k);
        
        % Piecewise constant plot using stairs
        stairs(x, control_plot(:,l), 'LineWidth', 1.5)
        
        ylim([0 1])
        xlim([0, (dn-1)*dx])
        xlabel('Position (mm)')
        ylabel('Control')
        title("weight = " + info{k,1} + ...
            ", c = " + info{k,2} + ...
            " \times c^*, k = " + info{k,3} + " \times k^*")
        
    end
    drawnow nocallbacks
    pause(0.001)
end




%% ANIMATE THE INPUTS BASED ON LUMINOSITY %%

figure(2)
clf
pause(.5)




for l=1:N

    for k=1:lidx

        subplot(length(idx),1,k)
        control_plot = controls{idx(k)}/control_norm(k);

        % Piecewise constant plot using stairs
        %
        imagesc(control_plot(:,l)',[0,1])
        colormap(gray); 
        axis image;
        
        xlabel('Position (pxl)')
        
        title("weight = " + info{k,1} + ...
            ", c = " + info{k,2} + ...
            " \times c^*, k = " + info{k,3} + " \times k^*.")
        
    end
    
    drawnow nocallbacks
    pause(0.001)

end





%% CREATE ANIMATED GIF %%

if createAnimatedGIF

    % --- GIF parameters ---
    
    datePresent = datetime('now','Format','yyyy-MM-dd_HH-mm');
    s = char(datePresent);
    
    gifname   = ['videoOutput/Optimal_inputs_1D_',s,'.gif'];
    delayTime = 0.05;   % seconds between frames
    firstFrame = true;
    
    if exist(gifname, 'file')
        delete(gifname)
    end
    
    
    % --- Figure setup ---
    
    fig = figure(3); 
    clf(fig,'reset') 
    set(gcf,'Position',[100 100 600 (length(idx)*100)])
    
    for l = 1:N
    
        for k = 1:length(idx)
            subplot(length(idx),1,k)
            control_plot = controls{idx(k)}/control_norm(k);
    
            imagesc(control_plot(:,l)', [0,1])
            colormap(gray)
            axis image
    
            xlabel('Position (pxl)')
            title("weight = " + info{k,1} + ...
                  ", c = " + info{k,2} + ...
                  " \times c^*, k = " + info{k,3} + " \times k^*")
        end
    
        drawnow
    
        % --- Capture one frame per l ---
        frame = getframe(gcf);
        im = frame2im(frame);
        [A,map] = rgb2ind(im,256);
    
        if firstFrame
            imwrite(A, map, gifname, 'gif', ...
                    'LoopCount', Inf, ...
                    'DelayTime', delayTime);
            firstFrame = false;
        else
            imwrite(A, map, gifname, 'gif', ...
                    'WriteMode', 'append', ...
                    'DelayTime', delayTime);
        end
    end
end


