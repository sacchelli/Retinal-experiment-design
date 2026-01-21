function [controls, inputParameters] = generate2DRandomControls(sdn, N, K, seed)
% GENERATE1DRANDOMCONTROLS Creates random control inputs for retinal experiments
%
% Inputs:
%   sdn             - Number of spatial discretization points in 1D
%   N               - Number of time steps
%   K               - Number of random controls to generate
%   seed            - Random seed for reproducibility
%
% Outputs:
%   controls        - Cell array {K x 1} of control matrices (dn x N)
%   inputParameters - Matrix (2 x K) with [NaN; NaN] for each control (no parameters)

% Handle optional arguments


minValue = 0;
maxValue = 1;


rng(seed);


% Pre-allocate
controls = cell(K, 1);
inputParameters = zeros(3, K);

for i = 1:K
    % Generate random control
    control_temp = minValue + (maxValue - minValue) * rand(sdn^2, N);
    
    controls{i} = control_temp;
    inputParameters(:, i) = [NaN; NaN; NaN];
end

end
