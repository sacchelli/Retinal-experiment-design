% The following script can be used to repoduce the numerical experiments
% present in the paper "...." by ....

% In the present script, we only select the optimal experiment on the basis
% of a nominal value, the maximum likelihood estimation in itself is not
% achieved here.

% This particular script focuses on the 1-dimensional model of the retina,
% which requires slightly different adjustments in the code.



addpath(genpath('./Functions'))
clc
timeTakenStart = tic;



%% FRONT MATTER %%


dn = 100; % Size of one layer

N = 100; % Number of measurements


% Pick the rng seed if necessary.

seed = 1;

rng(seed);


%%% Should we cutoff correlations?

% If cutoff is 'true', we do not use the full matrices for correlations, 
% using the fact that there is a forgetting effect in the
% dynamics.

cutoff = true;

% If cutoff is 'true', then depth determines the forgetting time
% beyond which we assume we can discard the correlations. Experimenting,
% it seems  depth * q. 20 is a good choice under current data. It is
% possible to experiment with the function 'plotDiagonalMax' to see the
% decrease in amplitude of the matrices along the diagonal. Naturally, this
% effect is highly dependent on Delta.

depth = 20;
% depth = 15;

%%% Should we add random controls in addition to waves?

addRandomInputs = true;

% How many?

nbrRandomInputs = 100; % must be >0 if true.



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

q = dn; % number of outputs


delta = retinalWidth/dn; % distance between cells


%%% Determination of size of the noise

% Notice that here, if we keep id as the noise matrix then only the
% relative sizes of the noises matter in the end.
% D is the dynamical noise and Sigmap is the output/measurement noise.


% D = 0.1 * eye(n); 
% Sigmap = 0.01 * eye(q); 

% D = 0.1 * eye(n); 
% Sigmap = 0.1 * eye(q); 

D = 0.1 * eye(n); 
Sigmap = 0.5 * eye(q); 


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

% 'absoluteLaplacian' connection graph:
%           Bi
%         / | \
%        1  2  1
%       /   |   \
%   A(i-1)  Ai   A(i+1)

% 'Laplacian' connection graph:
%           Bi
%         / | \
%        1  -2 1
%       /   |   \
%   A(i-1)  Ai   A(i+1)

% 'closestNeighbor' connection graph:
%           Bi
%         / | \
%        1  0  1
%       /   |   \
%   A(i-1)  Ai   A(i+1)


% Gamma = laplacian1DGraph(dn);
% Gamma = absoluteLaplacian1DGraph(dn);
% Gamma = laplacian1DGraph(dn);
 Gamma = closestNeighbor1DGraph(dn);


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

Acell{4} = [zz , -Gamma, zz ; zz, zz, zz ; zz, zz, zz];
Acell{5} = [zz , zz, zz ; Gamma, zz, zz ; zz, zz, zz];

Acell{6} = [zz , zz, zz ; zz, zz, zz ; GammaBG, zz, zz];
Acell{7} = [zz , zz, zz ; zz, zz, zz ; zz, GammaAG, zz];

Acell{8} = 1/tauB*Acell{1} + 1/tauA*Acell{2} + 1/tauG*Acell{3} ...
    + wm*Acell{4} + wp*Acell{5} + wBG*Acell{6} + wAG*Acell{7};

% Acell{8} contains the experimental value, the others serve as variations
% around that value.


Delta = max([tauA,tauB,tauG])/5; % Time between two measurements

Tf = N * Delta; % Total duration of the experiment

prec = 500;   % must be even



%%% Input and output matrices

B = sparse([eye(dn,m);zeros(2*dn,m)]);

C = sparse([zeros(q,2*dn),eye(q,dn)]);


%%% Contraction of A

thetaData = [1/tauB, 1/tauA, 1/tauG, wm, wp, wBG, wAG]; % Compilation of the experimental data

% We can choose the nominal vector theta, either as thetaData or as a
% variation of it. Again, we made the choice for this code that  the 
% values input in thetaNominal are variations around thetaData (already 
% contained in Acell{8}).

thetaNominal = 0.*thetaData;
%thetaNominal = (rand(1,p)-1/2)/2.*thetaData;


A = Acell{p+1};
for i = 1:p
    A = A + thetaNominal(i) * Acell{i};
end

maxReEigA = max(real(eigs(A)));



fprintf(['Data created, size of one layer: ', num2str(dn),', Max Re(λ(A)) = ', num2str(maxReEigA), '. \n'])

 
%% PRECOMPUTATION STEP %%

%%% Construction of sensitivities

% This construction follows the paper precisely.


fprintf(['\n','Dynamics and variance pre-computation... '])

tStepStart = tic;

F = buildF(A, Delta);
G = buildG(A, B, Delta, prec);
Sigma = buildSigma(A, D, Delta, prec);

% Use a cutoff accelerate the computations by killing off diagonal terms:

if cutoff
    SigmaBold = buildSigmaBoldCutoff(C, C/2, F, Sigma, Sigmap, N,depth);
    SigmaBoldinv = bandBlockDiagonalFilter(inv(full(SigmaBold)),depth,q,q); 
else
    SigmaBold = buildSigmaBold(C, C/2, F, Sigma, Sigmap, N);
    SigmaBoldinv = inv(SigmaBold);
end


%% PARTIAL DERIVATIVES PRE-COMPUTATIONS %%

tStepEnd = toc(tStepStart);

fprintf(['done (', num2str(tStepEnd,3) ,'s). \n'])

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

    Fb = buildF(Ab, Delta);
    
    Gb = buildG(Ab, Bb, Delta, prec);
    
    Sigmab = buildSigma(Ab, Db, Delta, prec);
    

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

%% COLLECTION OF INPUTS %%

%%% Test family: waves 

fprintf(['\n','Input collection... '])

tStepStart = tic;

cStar = retinalWidth/(100*Delta);
% cStar = retinalWidth/(Tf);

cMax = 2*cStar;   % max wave velocity
cMin = -2*cStar; % min wave velocity
% cMin = 0; % min wave velocity

kMax = 2*pi/retinalWidth*(dn-1)/2;   % max wave number
kMin = 2*pi/(retinalWidth);   % min wave number

cN = 31;
kN = 30;

[controls1DWave, inputParameters1Dwave] = generate1DPlaneWaveControls(dn, N, Delta, Tf, retinalWidth, cMin, cMax, cN, kMin, kMax, kN);

if addRandomInputs % produce random inputs if asked
    [controlsRandom, inputParametersRandom] = generate1DRandomControls(dn, N, nbrRandomInputs, seed);
end



% Concatenate inputs and constant input.

if addRandomInputs 
    counter = length(controls1DWave) + length(controlsRandom) + 1;
    inputParameters = [inputParameters1Dwave,inputParametersRandom,[NaN;0]];
    controls = cell(counter,1);
    controls(1:end-1) = [controls1DWave;controlsRandom];
else
    counter = length(controls1DWave) + 1;
    inputParameters = [inputParameters1Dwave,[NaN;0]];
    controls = cell(counter,1);
    controls(1:end-1) = controls1DWave;
end

% 

controls{counter} = ones(m, N); % Constant input at the end


K = counter;

tStepEnd = toc(tStepStart);

fprintf(['done (',num2str(tStepEnd,3) ,'s) \n'])

fprintf(['There are ', num2str(K),' inputs in this experiment'])



%% FISHER MATRICES CONSTRUCTION %%

% Please notice that the Fisher matrices are rescaled 

%%% Collection of Fisher matrices

fprintf(['\n','Fisher Matrices computations... \n'])

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

    M{i} = Mi/(N*q); % Rescale of final Fisher information matrix in order to account for influece of N and q
    
    if mod(i,floor(K/10))==0
        fprintf(['Progress:  %3.0f%% '],ceil(100*i/K))
        tStepEnd = toc(tStepStart);
        fprintf(['(',num2str(tStepEnd,3),' s).\n'])
        tStepStart = tic;
    end
end


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

Mtop = M{maxScoreIndex};


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

maxStepSearch = 1000; % Max number of loops
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


%% Animate the input

figure(1)
clf
pause(1)

% Define spatial domain

dx = retinalWidth/(dn-1);
x = 0:dx:(dn-1)*dx;

lidx = length(idx);
info = cell(lidx,3);
for k=1:length(idx)
    info{k,1} = num2str(wnew(idx(k)));
    info{k,2} = num2str(inputParameters(1,idx(k))/cStar);
    info{k,3} = num2str(inputParameters(2,idx(k))/kMax);
end

for l=1:N
    for k=1:length(idx)
        subplot(length(idx),1,k)
        control_plot = controls{idx(k)};
        
        % Piecewise constant plot using stairs
        stairs(x, control_plot(:,l), 'LineWidth', 1.5)
        
        ylim([0 2])
        xlim([0, (dn-1)*dx])
        xlabel('Position (mm)')
        ylabel('Control')
        title("weight = " + info{k,1} + ...
            ", c = " + info{k,2} + ...
            " \times c^*, k = " + info{k,3} + " \times k^*")
        
    end
    drawnow
    pause(0.001)
end




%% Animate based on luminosity

figure(2)
clf
pause(1)



lidx = length(idx);
info = cell(lidx,3);
for k=1:length(idx)
    info{k,1} = num2str(wnew(idx(k)));
    info{k,2} = num2str(inputParameters(1,idx(k))/cStar);
    info{k,3} = num2str(inputParameters(2,idx(k))/kMax);
end

for l=1:N
    for k=1:length(idx)
        subplot(length(idx),1,k)
        control_plot = controls{idx(k)};

        % Piecewise constant plot using stairs
        %
        imagesc(control_plot(:,l)',[0,2])
        colormap(gray); 
        axis image;
        
        xlabel('Position (pxl)')
        
        title("weight = " + info{k,1} + ...
            ", c = " + info{k,2} + ...
            " \times c^*, k = " + info{k,3} + " \times k^*")
        
    end
    drawnow
    pause(0.001)

end

%% Create animated gif %%

% --- GIF parameters ---

date_present = datetime('now','Format','yyyy-MM-dd_HH-mm');
s = char(date_present);

gifname   = ['video_output/Optimal_inputs_1D_',s,'.gif'];
delayTime = 0.05;   % seconds between frames
firstFrame = true;

if exist(gifname, 'file')
    delete(gifname)
end


% --- Figure setup ---

figure(3); 
clf(3,'reset') 
set(gcf,'Position',[100 100 600 (length(idx)*100)])

for l = 1:N

    clf(3)   % clear figure at each l

    for k = 1:length(idx)
        subplot(length(idx),1,k)
        control_plot = controls{idx(k)};

        imagesc(control_plot(:,l)', [0,2])
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
