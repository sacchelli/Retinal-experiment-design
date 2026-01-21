function M = gaussianPooling2DGraph(dn, delta, sigma)
% GAUSSIANPOOLING2DGRAPH Constructs a 2D Gaussian similarity (RBF) matrix
%
% Generates a dense matrix representing the Gaussian correlation between 
% all pairs of nodes in a 2D grid of size dn x dn.
%
% The formula for the interaction between node k (at coordinates xk, yk) 
% and node l (at coordinates xl, yl) is:
%   M(l,k) = exp( -((xk-xl)^2 + (yk-yl)^2) * delta^2 / (2 * sigma^2) )
%
% Inputs:
%   dn    - Number of points per dimension (total nodes = dn^2)
%   delta - Grid spacing (distance between adjacent nodes)
%   sigma - Standard deviation of the Gaussian kernel
%
% Output:
%   M     - Dense dn^2 x dn^2 Gaussian similarity matrix

% Precompute the base exponential factor for a distance of 1 unit
e1 = exp(-delta^2/(2*sigma^2));
M = zeros(dn^2, dn^2);

for k = 1:dn^2
    % Convert linear index k to 2D grid coordinates (ik, jk)
    % Note: mod(k,dn) results in 0 for the last element of a column
    ik = mod(k-1, dn) + 1; 
    jk = floor((k-1)/dn) + 1;
    
    for l = 1:dn^2
        % Convert linear index l to 2D grid coordinates (il, jl)
        il = mod(l-1, dn) + 1;
        jl = floor((l-1)/dn) + 1;
        
        % Euclidean distance squared: dist^2 = (ik-il)^2 + (jk-jl)^2
        % M(l,k) = (e1)^(dist^2)
        M(l,k) = e1^((ik-il)^2 + (jk-jl)^2);
    end
end


% %%% Potential short but less explainable alternative:
% M1D = gaussianPooling1DGraph(dn, delta, sigma);
% 
% % The 2D Gaussian kernel on a grid is the Kronecker product
% Mbis = kron(M1D, M1D);

end