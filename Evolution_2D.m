addpath(genpath('./Functions'))
clc
timeTakenStart = tic;

%% FRONT MATTER %%

sdn = 30; % Size of the edge of one layer
dn = sdn^2; % Size of one layer 

N = 100; % Number of measurements

%%% Should we cutoff correlations?
cutoff = true;
depth = 20;

%%% Should we add random controls in addition to waves?
addRandomInputs = false;
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

delta = retinalWidth/sdn; % distance between cells

%%% Determination of size of the noise
D = 0.1 * eye(n); 
Sigmap = 0.5 * eye(q); 

%%%% Model topology - CREATE ALL 4 TOPOLOGIES

fprintf('Building topologies...\n')

% Create all 4 topology matrices
Gamma1 = adjacencyLaplacian2DGraph(sdn);
Gamma2 = absoluteLaplacian2DGraph(sdn);
Gamma3 = laplacian2DGraph(sdn);
Gamma4 = closestNeighbor2DGraph(sdn);

topologyNames = {'Adjacency', 'Absolute', 'Laplacian', 'Closest'};

% B -> G and A -> G connections
GammaBG = gaussianPooling2DGraph(sdn,delta,sigmaBG);
GammaAG = gaussianPooling2DGraph(sdn,delta,sigmaAG);

%%% Construction of the matrix A for each topology
id = speye(dn);
zz = sparse(dn,dn);

% Store A matrices for all 4 topologies
F_all = cell(4,1);
G_all = cell(4,1);

thetaData = [1/tauB, 1/tauA, 1/tauG, wm, wp, wBG, wAG];
thetaNominal = 0.*thetaData;

Delta = max([tauA,tauB,tauG])/5; % Time between two measurements
Tf = N * Delta; % Total duration of the experiment
prec = 500;   % must be even

%%% Input and output matrices
B = sparse([eye(dn,m);zeros(2*dn,m)]);
C = sparse([zeros(q,2*dn),eye(q,dn)]);

% Build dynamics for all 4 topologies
Gammas = {Gamma1, Gamma2, Gamma3, Gamma4};

fprintf('Building dynamics matrices (this takes a moment)...\n')

for top = 1:4
    fprintf('  Topology %d/%d...', top, 4)
    tic_top = tic;
    
    Gamma = Gammas{top};
    
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
    
    A = Acell{p+1};
    for i = 1:p
        A = A + thetaNominal(i) * Acell{i};
    end
    
    F_all{top} = buildF(A, Delta);
    G_all{top} = buildG(A, B, Delta, prec);
    
    fprintf(' done (%.2fs)\n', toc(tic_top))
end

%% COLLECTION OF INPUTS %%

fprintf('\nGenerating inputs... ')
tStepStart = tic;

cStar = retinalWidth/(100*Delta);
cMax = 2*cStar;
cMin = 0;

kMax = 2*pi/retinalWidth*(sdn-1)/2;
kMin = 2*pi/(retinalWidth);

cN = 10;
kN = 10;
thetaN = 16;

[controls2DWave, inputParameters2Dwave] = generate2DPlaneWaveControls(sdn, N, Delta, retinalWidth, thetaN, cMin, cMax, cN, kMin, kMax, kN);

if addRandomInputs
    [controlsRandom, inputParametersRandom] = generate2DRandomControls(sdn, N, nbrRandomInputs, seed);
end

% Concatenate inputs
if addRandomInputs 
    counter = length(controls2DWave) + length(controlsRandom) + 1;
    inputParameters = [inputParameters2Dwave,inputParametersRandom,[NaN;NaN;0]];
    controls = cell(counter,1);
    controls(1:end-1) = [controls2DWave;controlsRandom];
else
    counter = length(controls2DWave) + 1;
    inputParameters = [inputParameters2Dwave,[NaN;NaN;0]];
    controls = cell(counter,1);
    controls(1:end-1) = controls2DWave;
end

controls{counter} = ones(m, N); % Constant input at the end

K = counter;

fprintf('done (%.2fs)\n', toc(tStepStart))
fprintf('Total inputs: %d\n', K)

%% SELECT INPUT AND SETUP FIGURE



% controlID = randi(K);
controlID = 1;
theta = 1.18*pi; % wave direction of propagation
c = 13*cStar; % wave velocity
kNumber = 5; % wave number

control_temp2D = zeros(sdn,sdn);
control_temp1D = zeros(sdn*sdn, N);
for it = 1:N
    t = (it-1) * Delta;
    for ix = 1:sdn
        x = (ix-1) * dx;  % Start at x=0
        for iy = 1:sdn
            y = (iy-1) * dx; % Start at y=0
            control_temp2D(iy,ix) = 1 + cos(kNumber*(cos(theta)*x + sin(theta)*y) - c*kNumber*t);
        end 
    end
    control_temp1D(:,it) = control_temp2D(:);
end



fprintf('\nUsing control ID: %d\n', controlID)
%u = controls{controlID};

u = control_temp1D;

fprintf(['Propagation angle: ',num2str(inputParameters(1,controlID)/pi),' pi\n'])

% Setup figure BEFORE simulation
fprintf('Setting up visualization...\n')
figure(1)
clf
set(gcf, 'Position', [50, 50, 1600, 900])

% Define spatial domain
dx = retinalWidth/(sdn-1);

% Layout: 4 rows (input, B, A, G) × 4 columns (4 topologies)
nRows = 4;
nCols = 4;

rowNames = {'Input', 'B', 'A', 'G'};

% Pre-create all subplots
ax_handles = gobjects(nRows, nCols);
im_handles = gobjects(nRows, nCols);

for row = 1:nRows
    for col = 1:nCols
        subplot_idx = (row-1)*nCols + col;
        ax_handles(row, col) = subplot(nRows, nCols, subplot_idx);
        
        % Create dummy image to get handle
        im_handles(row, col) = imagesc(zeros(sdn, sdn), 'Parent', ax_handles(row, col));
        axis(ax_handles(row, col), 'image')
        set(gca,'YDir','normal')
        % Set colormap
        if row == 1
            colormap(ax_handles(row, col), gray)
        else
            colormap(ax_handles(row, col), parula)
        end
        
        % Add title only to top row
        if row == 1
            title(ax_handles(row, col), topologyNames{col}, 'FontWeight', 'bold')
        end
        
        if col == 1
            ylabel(ax_handles(row, col), rowNames{row})
        else
            set(ax_handles(row, col), 'YTickLabel', [])
        end
        
        if row == 4
            xlabel(ax_handles(row, col), 'Position')
        else
            set(ax_handles(row, col), 'XTickLabel', [])
        end
    end
end

title_handle = sgtitle('Initializing...', 'FontSize', 14, 'FontWeight', 'bold');
drawnow

%% SIMULATE AND ANIMATE IN REAL-TIME

fprintf('Starting simulation and animation...\n')

% Initialize states for all topologies
x_all = zeros(n, 4);

for l = 1:N
    if mod(l, 10) == 0
        fprintf('  Step %d/%d\n', l, N)
    end
    
    ut = u(:,l);
    input_im = reshape(ut, sdn, sdn);
    
    % Update all topologies
    for top = 1:4
        x_all(:, top) = F_all{top} * x_all(:, top) + G_all{top} * ut;
        
        % Extract states
        stateB = reshape(x_all(1:dn, top), sdn, sdn);
        stateA = reshape(x_all(dn+1:2*dn, top), sdn, sdn);
        stateG = reshape(x_all(2*dn+1:3*dn, top), sdn, sdn);
        
        % Update images
        set(im_handles(1, top), 'CData', input_im);
        set(im_handles(2, top), 'CData', stateB);
        set(im_handles(3, top), 'CData', stateA);
        set(im_handles(4, top), 'CData', stateG);
    end
    
    % Update title
    set(title_handle, 'String', sprintf('Time step: %d/%d', l, N));
    pause(0.001)
    drawnow limitrate
end

fprintf(['Propagation angle: ',num2str(inputParameters(1,controlID)/pi),' pi'])
fprintf('\nAnimation complete!\n') 
fprintf('Total time: %.2fs\n', toc(timeTakenStart))