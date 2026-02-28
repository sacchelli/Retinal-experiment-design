function M = gaussianPooling1DGraph(N, delta, sigma)
% GAUSSIANPOOLING1DGRAPH Constructs a Gaussian radial basis function (RBF) matrix
%
% Generates an N x N symmetric matrix where each entry M(i,j) represents the 
% Gaussian correlation between nodes i and j based on their distance. 
%
% The formula used for each entry is:
%   M(i,j) = exp( -(|i-j| * delta)^2 / (2 * sigma^2) )
%
% Inputs:
%   N     - Number of points in the 1D graph
%   delta - Spatial distance between adjacent nodes (grid spacing)
%   sigma - Standard deviation (controls the "width" or "reach" of the pooling)
%
% Output:
%   M     - Dense N x N Gaussian similarity matrix

val = exp(-delta^2/(2*sigma^2));
M = eye(N); 
for i = 1:N-1
    % Compute the scalar value for this specific diagonal (lag i)
    dist_sq = i^2;
    M_val = val^(dist_sq);
    % Fill the diagonals with that scalar
    M = M + diag(M_val*ones(N-i,1), i) + diag(M_val*ones(N-i,1), -i);
end
M = M/(sqrt(2*pi)*sigma);
end