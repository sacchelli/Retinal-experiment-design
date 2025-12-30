function [controls, inputParameters] = generate1DRandomControls(dn, N, K, seed)
% GENERATE1DRANDOMCONTROLS Creates random control inputs for retinal experiments
%
% Inputs:
%   dn              - Number of spatial discretization points
%   N               - Number of time steps
%   K               - Number of random controls to generate
%   seed            - Random seed for reproducibility
%
% Outputs:
%   controls        - Cell array {K x 1} of control matrices (dn x N)
%   inputParameters - Matrix (2 x K) with [NaN; NaN] for each control (no parameters)

% Handle optional arguments


minValue = 0;
maxValue = 2;
spatialSmoothing = 0;
temporalSmoothing = 0;

rng(seed);


% Pre-allocate
controls = cell(K, 1);
inputParameters = zeros(2, K);

for i = 1:K
    % Generate random control
    control_temp = minValue + (maxValue - minValue) * rand(dn, N);
    
    % Apply spatial smoothing if requested
    if spatialSmoothing > 0
        for t = 1:N
            control_temp(:, t) = gaussianSmooth1D(control_temp(:, t), spatialSmoothing);
        end
    end
    
    % Apply temporal smoothing if requested
    if temporalSmoothing > 0
        for x = 1:dn
            control_temp(x, :) = gaussianSmooth1D(control_temp(x, :), temporalSmoothing);
        end
    end
    
    % Ensure values stay in bounds after smoothing
    control_temp = max(minValue, min(maxValue, control_temp));
    
    controls{i} = control_temp;
    inputParameters(:, i) = [NaN; NaN];
end

end

% Helper function to get field with default value
function val = getfield_default(s, field, default)
    if isfield(s, field)
        val = s.(field);
    else
        val = default;
    end
end

% Helper function for 1D Gaussian smoothing
function smoothed = gaussianSmooth1D(data, sigma)
    if sigma <= 0
        smoothed = data;
        return;
    end
    
    % Create Gaussian kernel
    kernelSize = ceil(6 * sigma);
    if mod(kernelSize, 2) == 0
        kernelSize = kernelSize + 1;
    end
    x = -(kernelSize-1)/2 : (kernelSize-1)/2;
    kernel = exp(-x.^2 / (2*sigma^2));
    kernel = kernel / sum(kernel);
    
    % Apply convolution with symmetric boundary conditions
    smoothed = conv(data(:), kernel, 'same');
    smoothed = reshape(smoothed, size(data));
end