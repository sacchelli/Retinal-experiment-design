% The following script can be used to reproduce the numerical experiments
% present in the paper "...." by ....

% This particular script focuses on the 1-dimensional model of the retina,
close("all")
%clear variables

% Colors
Goldenrod=[1 0.9 0.15];
red=[0.9,0.0,0.0];
yellow="#FAB84D";
blue="#3DA3FF";
green="#30E172";
purple="#8D4AC8";
grey1="#F2F2F2";
grey2="#C8C8C8";
grey1="#848484";
grey4="#585858";

addpath(genpath('./Functions'))
clc
timeTakenStart = tic;

%% FRONT MATTER %%
%---------------------------------------------------------------------------------
% number of parameters: set p to 6 or 7 (different identifiability properties)
p = 7; 
% regularisation parameter for Fisher information matrix (e.g., 1e-6)
regul=0; 

% Comparison of computational times of Fisher information matrices?
%Compare_FIM_times=false;
Compare_FIM_times=true;

% Generate data with a deterministic model (without process noise and x0=0) 
%   or a stochastic one (with process noise and x0 random)?
    data_model='stochastic'
    %data_model='deterministic'

% Estimate parameters assuming a deterministic model (without process noise and x0=0) 
%   or a stochastic one (with process noise and x0 random)?
    estimation_model='stochastic'
    %estimation_model='deterministic'

% Generate design assuming a deterministic model (without process noise and x0=0) 
%   or a stochastic one (with process noise and x0 random)?
    design_model='stochastic'
    %design_model='deterministic'
% Optimality criterion
    crit = 'Dopt'
         KK=NaN; cc=NaN;
    %----------------------------------------------------------------------
    % crit = 'Lopt'; % two sub-cases
        %------------------------------------------------------------------
    %     crit_Lopt = 'Aopt'; 
    %        KK=eye(p); % A-optimal design
    %        cc=NaN;
    %        relative_prec=true; % true for relative precision (KK=diag(theta.^(-1)))
    %        %relative_prec=false; 
        %------------------------------------------------------------------
        % crit_Lopt = 'copt'; 
        %     KK=zeros(p,1); i=1; KK(i)=1; cc=KK; % c-optimal design for parameter nb. i
        %     %KK=[0.0001   -0.0001   -0.0000   -0.4810    0.7075    0.0000   -0.5177]';
        %      cc=KK;
        %------------------------------------------------------------------
required_eff=0.999; % required efficiency in the optimal design algorithm

% Generate an optimal design with thetaDesign=thetaData ?
    true_theta=true; 
% otherwise, thetaDesign=false_theta;  
    false_theta=0.01*ones(1,p);

% Remove inessential matrices ?   
    remove_inessential=true
    %remove_inessential=false
% number of experiments (= total number of inputs used)    
Nexp=10; 

%dn = 100; % Size of one layer
dn=20;
%dn=4;

N = 100; % Number of measurements
%N=50;
%N=2;

%---------------------------------------------------------------------------------


% Pick the rng seed if necessary.
seed = 5;
rng(seed);

%%% Should we cutoff correlations?
% If cutoff is 'true', we do not use the full matrices for correlations, 
% using the fact that there is a forgetting effect in the
% dynamics.

%cutoff = true;
cutoff = false;

% If cutoff is 'true', then depth determines the forgetting time
% beyond which we assume we can discard the correlations. Experimenting,
% it seems  depth 20 is a good choice under current data. It is
% possible to experiment with the function 'plotDiagonalMax' to see the
% decrease in amplitude of the matrices along the diagonal. Naturally, this
% effect is highly dependent on Delta.

depth = 20; % Seems to be a sweet spot under current values, 
            % increase Delta to decrease the depth.

%%% Should we add random controls in addition to waves?
addRandomInputs = false;

% How many?
nbrRandomInputs = 300; % must be >0 if true.

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
% size of the system
    n = 3 * dn; 
% number of inputs
    m = dn; 
% number of outputs
%   if q<dn, the observed responses are spread among the last dn components of the state
    %q = dn; 
    q=4; 
    positions = [floor(dn/5),floor(2*dn/5),floor(3*dn/5),floor(4*dn/5)]

delta = retinalWidth/dn; % distance between cells

%%% Determination of size of the noise
% Notice that here, if we keep id as the noise matrix then only the
% relative sizes of the noises matter in the end ---> No: the covariance
% appears differently in Q and the part that depends on inputs
% D is the dynamical noise and Sigmap is the output/measurement noise.

sigma_measurement=1;
beta=1;
%
sigma_process=beta*sigma_measurement;
D = sigma_process * eye(n); 
%
sig2=sigma_measurement^2;
Sigmap = sig2 * eye(q); 

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

% 'closestNeighbor' connection graph:
%           Bi
%         / | \
%        1  0  1
%       /   |   \
%   A(i-1)  Ai   A(i+1)

% 'Laplacian' connection graph (unrealistic for the retina):
%           Bi
%         / | \
%        1  -2 1
%       /   |   \
%   A(i-1)  Ai   A(i+1)

% Gamma = laplacian1DGraph(dn);
% Gamma = absoluteLaplacian1DGraph(dn);
Gamma = closestNeighbor1DGraph(dn);
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
Acell{4} = [zz , -Gamma, zz ; zz, zz, zz ; zz, zz, zz];
Acell{5} = [zz , zz, zz ; Gamma, zz, zz ; zz, zz, zz];

if p==7
    Acell{6} = [zz , zz, zz ; zz, zz, zz ; GammaBG, zz, zz];
    Acell{7} = [zz , zz, zz ; zz, zz, zz ; zz, -GammaAG, zz];
    
    Acell{8} = 1/tauB*Acell{1} + 1/tauA*Acell{2} + 1/tauG*Acell{3} ...
        + wm*Acell{4} + wp*Acell{5} + wBG*Acell{6} + wAG*Acell{7};
    % Acell{8} contains the experimental value, the others serve as variations
    % around that value.
elseif p==6
    Acell{6} = [zz , zz, zz ; zz, zz, zz ; GammaBG, -GammaAG, zz];
    
    Acell{7} = 1/tauB*Acell{1} + 1/tauA*Acell{2} + 1/tauG*Acell{3} ...
        + wm*Acell{4} + wp*Acell{5} + wBG*Acell{6};    
    % Acell{7} contains the experimental value, the others serve as variations
    % around that value.
end

%Delta = max([tauA,tauB,tauG])/2; % Time between two measurements
Delta = 8.25;

Tf = N * Delta; % Total duration of the experiment
%prec = 500;   % must be even

%%% Input and output matrices
B = sparse([eye(dn,m);zeros(2*dn,m)]);
% C = sparse([zeros(q,2*dn),eye(q,dn)]);
% % q points with maximum covering in [-1,1]
%     xCR=-1+1/q:(2/q):1-1/q;
%     % now in [0,1]
%     xCR=(xCR+1)/2;
%     % now integers spread in 1:dn
%     nCR=round((dn-1)*xCR+1);
%     C=zeros(q,3*dn);
%     for i=1:q
%         if nCR(i) <= dn % Vérifie que l'indice est valide
%             C(i, 2*dn+nCR(i)) = 1;
%         end
%     end

C = sparse([zeros(dn,2*dn),eye(dn,dn)]);
C = C(positions,:);   % defined before (where q is set) 
%%% Contraction of A
if p==7
    thetaData = [1/tauB, 1/tauA, 1/tauG, wm, wp, wBG, wAG]; % Compilation of the experimental data
elseif p==6
    thetaData = [1/tauB, 1/tauA, 1/tauG, wm, wp, wBG];
end

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

fprintf(['Data created, size of one layer: ', num2str(dn),', Max Re(λ(A)) = ', num2str(maxReEigA), '. \n'])
 

%% COLLECTION OF INPUTS %%

[controls, counter, cStar, kMax, inputParameters] = inputs_1D(m, dn, N, Delta, Tf, retinalWidth,...
                                addRandomInputs, nbrRandomInputs, seed );
K = counter;

%% FISHER MATRICES CONSTRUCTION %%
% 
% Please notice that the Fisher information matrices are rescaled 
%

if true_theta
    thetaDesign=thetaData;
else
    thetaDesign=false_theta;
end    

if Compare_FIM_times
    % Comparison of computational times (for 10 inputs only): 
    %   --> direct calculation/recursive calculation (Kalman type)
    
    t1=tic;
    M = MFisher(design_model,controls(1:10),thetaDesign, sig2,Delta, q, Acell, B, C,D, N,cutoff,depth,regul);
    t_Fisher_expand=toc(t1);
    M1=M;
    
    t1=tic;
    M = MFisher_recursive(design_model,controls(1:10),thetaDesign, sig2,Delta, q, Acell, B, C,D, N);
    t_Fisher_recursive=toc(t1);
    M2=M;
    
    time_ratio=t_Fisher_recursive/t_Fisher_expand;
    fprintf(['Ratio of computational times direct/recursive: ', num2str(time_ratio), '. \n'])
end

(M1{1}-M2{1})./M1{1}

return

% Calculation of all Fisher information matrices for parameters thetaData (by the fastest method):
t1=tic;
M = MFisher(design_model,controls,thetaData, sig2,Delta, q, Acell, B, C,D, N,cutoff,depth,regul);
t_Fisher=toc(t1);
fprintf(['Computational time of the ', num2str(K), ' matrices: ', num2str(t_Fisher,3), 's. \n'])

MDesign=M; % Matrices used to construct optimal designs
if thetaDesign~=thetaData
    % We construct now designs for parameters thetaDesign 
    % --> compute all Fisher information matrices for thetaDesign
    MDesign = MFisher(design_model,controls,thetaDesign, sig2,Delta, q, Acell, B, C,D, N,cutoff,depth,regul);
end


%% OPTIMAL INPUT SELECTION (for MDesign) %%
tSearchStart = tic;

if strcmp(crit,'Lopt') && strcmp(crit_Lopt,'Aopt') && relative_prec
    KK=diag(thetaDesign.^(-1));
end    
[Mopt,wopt,indices] = optimal_FIM(MDesign,[],crit,KK,20,required_eff,remove_inessential,100);
[nn,idx] = efficient_rounding(wopt,Nexp)

wnew=wopt(idx);
idx=indices(idx) % indices corresponding to the original indexation of matrices in M

%--------------------------------------------------------------------------    

tSearchTaken = toc(tSearchStart)

timeTakenCheck = toc(timeTakenStart);

fprintf(['\n','Computation time for optimal design: ', num2str(tSearchTaken,3) ,' s\n'])
fprintf(['\n','Total computation time up to now: ', num2str(timeTakenCheck,3) ,' s\n'])
fprintf([' \n'])
 
%% Finding the best vertex before initialization

maxScore = -Inf;
maxScoreIndex = 0;

for k = 1:K
    if maxScore < logdet(M{k})
        maxScore = logdet(M{k});
        maxScoreIndex = k;
    end
end

Mtop = M{maxScoreIndex};

% ANIMATE THE INPUTS %%

figure(1)
clf
pause(.5)

% Define spatial domain

dx = retinalWidth/(dn-1);
x = 0:dx:(dn-1)*dx;

lidx = length(idx);
info = cell(lidx,3);
for k=1:length(idx)
    %info{k,1} = num2str(wnew(idx(k)));
    info{k,1} = num2str(wnew(k));
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
    drawnow nocallbacks
    pause(0.001)
end

%return


%% Confirmation of the identifiability problem for the model with 7 parameters
if p==7 % (identifiability problem)

    DETM=NaN(1,K);
    for k=1:K
        DETM(k)=det(M{k});
    end
    idetneg=find(DETM<0);
    Ldetneg=length(idetneg)
    [detmin,imin]=min(DETM); detmin
    [V,Lambda]=eig(M{imin});
    [Lambdamin,iimin]=min(diag(Lambda)); 
    ee=diag(Lambda);
    ee(ee<1e-3)'
    V(:,iimin)'
    
    if Ldetneg>0 
        %for i=1:Ldetneg
        for i=1:2
            [V,Lambda]=eig(M{idetneg(i)});
            [Lambdamin,iimin]=min(diag(Lambda)); 
            ee=diag(Lambda);
            ee(ee<1e-3)'
            vv=V(:,iimin)'
            [vv(4)/vv(5) vv(7)/vv(5)]
            [-thetaDesign(4)/thetaDesign(5) -thetaDesign(7)/thetaDesign(5)] 
            pause(0.1)
        end
    end
end

%% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% Generate data for a design defined by idx (indices of inputs in {1,...,K=counter}
% and their associated weights

% Optimal inputs
% theoretical results (Fisher information matrix)
%[Cov_theoretical,sqrt_eig_Cov_theoretical] = Stat_THETA_theoretical(M,idx,nn,N,q);
Cov_theoretical = cov_THETA_theoretical(M,idx,nn,N,q);

if crit == 'Dopt'
    Cov_theoretical_Dopt=Cov_theoretical;
    MDopt=Mopt;
    Eff_Dopt=exp((1/p)*(-log(det(N*q*Nexp*Cov_theoretical)) - log(det(det(MDopt))) ));
elseif crit == 'Lopt'
    if crit_Lopt == 'Aopt'
        Cov_theoretical_Aopt=Cov_theoretical;
        MAopt=Mopt;
        Eff_Aopt=trace(inv(MAopt))/trace(N*q*Nexp*Cov_theoretical);
    elseif crit_Lopt == 'copt'
        Cov_theoretical_copt=Cov_theoretical;
        Mcopt=Mopt;
        Eff_copt=((cc'/Mopt)*cc)/(N*q*Nexp*cc'*Cov_theoretical*cc);
    end
end    
% Set variables as NaN when then don't exist
TestD=exist( 'Cov_theoretical_Dopt', 'var' );
if TestD == 0
    Cov_theoretical_Dopt=NaN; %NaN(p,p);
    Eff_Dopt=NaN;
end 
TestA=exist( 'Cov_theoretical_Aopt', 'var' );
if TestA == 0
    Cov_theoretical_Aopt=NaN; %NaN(p,p);
    Eff_Aopt=NaN;
end
Testc=exist( 'Cov_theoretical_copt', 'var' );
if Testc == 0
    Cov_theoretical_copt=NaN; %NaN(p,p);
    Eff_copt=NaN;
end
eig_Dopt=eig(Cov_theoretical_Dopt)';
eig_Aopt=eig(Cov_theoretical_Aopt)';
eig_copt=eig(Cov_theoretical_copt)';
diag_Dopt=diag(Cov_theoretical_Dopt)';
diag_Aopt=diag(Cov_theoretical_Aopt)';
diag_copt=diag(Cov_theoretical_copt)';

% Random inputs
idx_rand=counter-1-nbrRandomInputs:1:counter-1-nbrRandomInputs+Nexp; % indices for random inputs
nn_rand=ones(1,Nexp); % each one used one time only
%[Cov_theoretical_rand,sqrt_eig_Cov_theoretical_rand] = Stat_THETA_theoretical(M,idx_rand,nn_rand,N,q);
Cov_theoretical_rand = cov_THETA_theoretical(M,idx_rand,nn_rand,N,q);
eig_rand=eig(Cov_theoretical_rand)';
diag_rand=diag(Cov_theoretical_rand)';

efficiencies=[Eff_Dopt Eff_Aopt Eff_copt] % efficiencies of rounded designs wrt true optimal designs (which are 

% Summary:  rows    = Dopt, Aopt, copt, random
%           columns = Eff_D Eff_A Eff_A (with respect to rounded optimal designs)
Summary=[exp((log(det(Cov_theoretical_Dopt))-log(det(Cov_theoretical_Dopt)))/p) trace(Cov_theoretical_Aopt)/trace(Cov_theoretical_Dopt) trace(cc'*Cov_theoretical_copt*cc)/trace(cc'*Cov_theoretical_Dopt*cc)
         exp((log(det(Cov_theoretical_Dopt))-log(det(Cov_theoretical_Aopt)))/p) trace(Cov_theoretical_Aopt)/trace(Cov_theoretical_Aopt) trace(cc'*Cov_theoretical_copt*cc)/trace(cc'*Cov_theoretical_Aopt*cc) 
         exp((log(det(Cov_theoretical_Dopt))-log(det(Cov_theoretical_copt)))/p) trace(Cov_theoretical_Aopt)/trace(Cov_theoretical_copt) trace(cc'*Cov_theoretical_copt*cc)/trace(cc'*Cov_theoretical_copt*cc) 
         exp((log(det(Cov_theoretical_Dopt))-log(det(Cov_theoretical_rand)))/p) trace(Cov_theoretical_Aopt)/trace(Cov_theoretical_rand) trace(cc'*Cov_theoretical_copt*cc)/trace(cc'*Cov_theoretical_rand*cc)]

%%
figure(11)
pp1D=plot(1:p,log10(eig_Dopt),'-','LineWidth',2,'MarkerSize',10); hold on
pp1A=plot(1:p,log10(eig_Aopt),'-','LineWidth',2,'MarkerSize',10);
pp1c=plot(1:p,log10(eig_copt),'-','LineWidth',2,'MarkerSize',10);
pp1rand=plot(1:p,log10(eig_rand),'-','LineWidth',2,'MarkerSize',10);
set(pp1D,'Color',red,'LineStyle','--','Marker','p','MarkerFaceColor',red)
set(pp1A,'Color',blue,'LineStyle','--','Marker','o','MarkerFaceColor',blue)
set(pp1c,'Color',yellow,'LineStyle','--','Marker','d','MarkerFaceColor',yellow)
set(pp1rand,'Color','k','LineStyle','--','Marker','s','MarkerFaceColor','k')
set(gca,'FontSize',20)
tightfig(gcf);
axis tight
grid on
hold off

figure(12)
pp2D=plot(1:p,log10(diag_Dopt),'-','LineWidth',2,'MarkerSize',10); hold on
pp2A=plot(1:p,log10(diag_Aopt),'-','LineWidth',2,'MarkerSize',10);
pp2c=plot(1:p,log10(diag_copt),'-','LineWidth',2,'MarkerSize',10);
pp2rand=plot(1:p,log10(diag_rand),'-','LineWidth',2,'MarkerSize',10);
set(pp2D,'Color',red,'LineStyle','--','Marker','p','MarkerFaceColor',red)
set(pp2A,'Color',blue,'LineStyle','--','Marker','o','MarkerFaceColor',blue)
set(pp2c,'Color',yellow,'LineStyle','--','Marker','d','MarkerFaceColor',yellow)
set(pp2rand,'Color','k','LineStyle','--','Marker','s','MarkerFaceColor','k')
set(gca,'FontSize',20)
axis tight
grid on
hold off
tightfig(gcf);

return

%%

% -----------------------------------------------------------------------
% Generate only one data set (for thetaData) and estimate, just to check...
Y=data_deterministic_1D_rec(controls,idx,nn,sig2,thetaData, Delta, Acell, B, C, N);
    % Sufficient statistics: 
    [YbarD,CovYD] = average_data(Y,nn);

Y = data_stochastic_1D_rec(controls,idx,nn,sig2,thetaData, Delta, Acell, B, C, D, N);

    % Sufficient statistics: 
    [YbarS,CovYS] = average_data(Y,nn);
if strcmp(data_model,'deterministic')
    Ybar=YbarD; CovY=CovYD;
elseif strcmp(data_model,'stochastic')
    Ybar=YbarS; CovY=CovYS;
end

ntrials=1;
t1=tic;
for ii=1:ntrials
    JML1=JML_stochastic_direct(Ybar,CovY,controls,idx,nn,sig2,thetaData, Delta, Acell, B, C, D,cutoff, depth);
end
toc(t1)
t2=tic;
for ii=1:ntrials
    JML=JML_stochastic(Y,controls,idx,nn,sig2,thetaData, Delta, Acell, B, C, D);
end
toc(t2)

[JML1 JML]

return
%%

% Model response (deterministic 1D model) at parameters theta:
if strcmp(estimation_model,'deterministic') 
    % Weighted LS  
    J_ML_theta_1D=@(theta) JML_deterministic(Ybar,controls,idx,nn,sig2,theta, Delta, Acell, B, C);
        JML_0=JML_deterministic(Ybar,controls,idx,nn,sig2,thetaData, Delta, Acell, B, C)
    J_ML_grad_theta_1D=@(theta) JML_grad_deterministic(Ybar,controls,idx,nn,sig2,theta, Delta, Acell, B, C); 
elseif strcmp(estimation_model,'stochastic') 
    % ML
    J_ML_theta_1D=@(theta) JML_stochastic(Y,controls,idx,nn,sig2,theta, Delta, Acell, B, C, D);
        JML_0=JML_stochastic(Y,controls,idx,nn,sig2,thetaData, Delta, Acell, B, C, D)
end

%
options = optimoptions('fminunc', 'MaxFunctionEvaluations', 1e4, 'MaxIterations', 1e4, ...
    'Algorithm','quasi-newton', ...
    'SpecifyObjectiveGradient', false);
t1=tic;
[thetaopt_a,fval_a] = fminunc(J_ML_theta_1D,thetaData,options);         % optimisation initialised at thetaData
t1p=toc(t1);
% options = optimoptions('fminunc', 'MaxFunctionEvaluations', 1e4, 'MaxIterations', 1e4, ...
%     'Algorithm','quasi-newton', ...
%     'SpecifyObjectiveGradient', false); 
% This one uses gradient (only for the deterministic model!) and is much slower :
% options = optimoptions('fminunc', 'MaxFunctionEvaluations', 1e4, 'MaxIterations', 1e4, ...
%     'Algorithm', 'quasi-newton', ...
%     'SpecifyObjectiveGradient', true);
t2=tic;
%[thetaopt_b,fval_b] = fminunc(J_ML_grad_theta_1D,thetaData,options);    % optimisation initialised at thetaData
[thetaopt_b,fval_b] = fminunc(J_ML_theta_1D,thetaData,options);    % optimisation initialised at thetaData
t2p=toc(t2);
t3=tic;
%[thetaopt_c,fval_c] = fminunc(J_ML_grad_theta_1D,0.01*ones(1,p),options);    % optimisation initialised elsewhere (0.01*ones(1,p))
[thetaopt_c,fval_c] = fminunc(J_ML_theta_1D,0.01*ones(1,p),options);    % optimisation initialised elsewhere (0.01*ones(1,p))
t3p=toc(t3);
[thetaopt_a
 thetaopt_b
 thetaopt_c   
 thetaData]
[fval_a fval_b fval_c JML_0]
[t1p t2p t3p]

%% now repeat Nrepeat times...    
Nrepeat=0;
%Nrepeat=100;
THETA=NaN(Nrepeat,p);
timebar= waitbar(0,'repeated simulations & estimations...'); 
for i=1:Nrepeat
    waitbar(i/Nrepeat,timebar);
    if strcmp(data_model,'deterministic')
        Y=data_deterministic_1D_rec(controls,idx,nn,sig2,thetaData, Delta, Acell, B, C, N);
            % Sufficient statistics: 
            [Ybar,CovY] = average_data(Y,nn);
        if strcmp(estimation_model,'deterministic') 
            % Criterion for the estimation of theta: weighted LS 
            J_ML_theta_1D=@(theta) JML_deterministic(Ybar,controls,idx,nn,sig2,theta, Delta, Acell, B, C);        
        elseif strcmp(estimation_model,'stochastic')
            % Criterion for the estimation of theta: ML 
            J_ML_theta_1D=@(theta) JML_stochastic(Y,controls,idx,nn,sig2,theta, Delta, Acell, B, C, D);
        end
    elseif strcmp(data_model,'stochastic')
        Y = data_stochastic_1D_rec(controls,idx,nn,sig2,thetaData, Delta, Acell, B, C, D, N);
            [Ybar,CovY] = average_data(Y,nn);
        if strcmp(estimation_model,'deterministic') 
            % Criterion for the estimation of theta: weighted LS 
            J_ML_theta_1D=@(theta) JML_deterministic(Ybar,controls,idx,nn,sig2,theta, Delta, Acell, B, C);
        elseif strcmp(estimation_model,'stochastic') 
            % Criterion for the estimation of theta: ML 
            J_ML_theta_1D=@(theta) JML_stochastic(Y,controls,idx,nn,sig2,theta, Delta, Acell, B, C, D);
        end        
    end
    options = optimoptions('fminunc', 'MaxFunctionEvaluations', 1e4, 'MaxIterations', 1e4, ...
    'Algorithm','quasi-newton', ...
    'Display', 'none', ...
    'SpecifyObjectiveGradient', false);

    evalc('[thetaopt,fval_b] = fminunc(J_ML_theta_1D,0.01*ones(1,p))');
    %evalc('[thetaopt,fval_b] = fminunc(J_ML_theta_1D,thetaData)');
    THETA(i,:)=thetaopt;   
end    
close(timebar)

%--------------------------------------------------------------------------------------
    % Results summary
    format shortG
    [mean_theta,mean_absbias_theta,Cov_theta,sqrt_eig_Cov] = Stat_THETA(THETA,thetaData)
    format short
%--------------------------------------------------------------------------------------

return

%% ANIMATE THE INPUTS BASED ON LUMINOSITY %%

figure(2)
clf
pause(.5)



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
        imagesc(control_plot(:,l)',[0,1])
        colormap(gray); 
        axis image;
        
        xlabel('Position (pxl)')
        
        title("weight = " + info{k,1} + ...
            ", c = " + info{k,2} + ...
            " \times c^*, k = " + info{k,3} + " \times k^*")
        
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
            control_plot = controls{idx(k)};
    
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

