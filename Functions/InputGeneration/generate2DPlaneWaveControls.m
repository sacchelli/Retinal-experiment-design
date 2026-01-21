function [controls, inputParameters] = generate2DPlaneWaveControls(sdn, N, Delta, retinalWidth, thetaN, cMin, cMax, cN, kMin, kMax, kN)
% GENERATE1DPLANEWAVECONTROLS Creates 1D traveling plane wave controls for retinal experiments
%
% Inputs:
%   sdn             - Number of spatial discretization points in one
%                       dimension
%   N               - Number of time steps
%   Delta           - Time step size (ms)
%   Tf              - Total experiment time (ms)
%   retinalWidth    - Width of retina (mm)
%   thetaN          - Number of orientation values (rad)
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

dx = retinalWidth/(sdn-1);

% Generate velocity and wave number grids

thetaValues = linspace(0,2*pi*(1-1/thetaN),thetaN);
cValues = linspace(cMin, cMax, cN);
kValues = linspace(kMin, kMax, kN);

% Calculate total number of controls: cN*kN 
K = thetaN * cN * kN;

% Pre-allocate
controls = cell(K, 1);
inputParameters = zeros(3, K);

counter = 0;

% Generate traveling wave controls
for itheta = 1:thetaN
    for ic = 1:cN
        for ik = 1:kN
            counter = counter + 1;
            
            % Parameters
            theta = thetaValues(itheta); % wave direction of propagation
            c = cValues(ic); % wave velocity
            kNumber = kValues(ik); % wave number
            
            control_temp2D = zeros(sdn,sdn);
            control_temp1D = zeros(sdn*sdn, N);
            for it = 1:N
                t = (it-1) * Delta;
                for ix = 1:sdn
                    x = (ix-1) * dx;  % Start at x=0
                    for iy = 1:sdn
                        y = (iy-1) * dx; % Start at y=0
                        control_temp2D(iy, ix) = 0.5*(1 + cos(kNumber*(cos(theta)*x + sin(theta)*y) - c*kNumber*t)); % Careful about the position (iy,ix) to be consistent with the column-major ordering
                    end 
                end
                control_temp1D(:,it) = control_temp2D(:);
            end
            
            controls{counter} = control_temp1D;
            inputParameters(:, counter) = [theta; c; kNumber];
        end
    end
end
end