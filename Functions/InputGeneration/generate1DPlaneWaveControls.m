function [controls, inputParameters] = generate1DPlaneWaveControls(dn, N, Delta, Tf, retinalWidth, cMin, cMax, cN, kMin, kMax, kN)
% GENERATE1DPLANEWAVECONTROLS Creates 1D traveling plane wave controls for retinal experiments
%
% Inputs:
%   dn              - Number of spatial discretization points
%   N               - Number of time steps
%   Delta           - Time step size (ms)
%   Tf              - Total experiment time (ms)
%   retinalWidth    - Width of retina (mm)
%   cMin            - Minimum wave velocity (mm/ms)
%   cMax            - Maximum wave velocity (mm/ms)
%   cN              - Number of velocity values
%   kMin            - Minimum wave number (rad/mm)
%   kMax            - Maximum wave number (rad/mm)
%   kN              - Number of wave numbers
%
% Outputs:
%   controls        - Cell array {K x 1} of control matrices (dn x N)
%   inputParameters - Matrix (2 x K) with [velocity; wavenumber] for each control



% Setup spatial discretization

dx = retinalWidth/(dn-1);

% Generate velocity and wave number grids
cValues = linspace(cMin, cMax, cN);
kValues = linspace(kMin, kMax, kN);

% Calculate total number of controls: cN*kN 
K = cN * kN;

% Pre-allocate
controls = cell(K, 1);
inputParameters = zeros(2, K);

counter = 0;

% Generate traveling wave controls
for ic = 1:cN
    for ik = 1:kN
        counter = counter + 1;
        
        % Parameters
        c = cValues(ic); %wave velocity
        kNumber = kValues(ik); % wave number
        
        control_temp = zeros(dn, N);
        for it = 1:N
            t = (it-1) * Delta;
            for ix = 1:dn
                x = (ix-1) * dx;  % Start at x=0
                control_temp(ix, it) = 0.5*(1 + cos(kNumber*x - c*kNumber*t));
            end
        end
        
        controls{counter} = control_temp;
        inputParameters(:, counter) = [c; kNumber];
    end
end

end