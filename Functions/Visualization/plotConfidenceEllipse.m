function plotConfidenceEllipse(Cov2, mu, level, color)
% plotConfidenceEllipse(Cov2, mu, level, color)
% Plots a 2D confidence ellipse for a 2x2 covariance matrix Cov2,
% centred at mu = [mu1, mu2], at confidence level (e.g. 0.95).

% Chi-squared threshold for 2 degrees of freedom
chiSq = -2 * log(1 - level);   % = chi2inv(level, 2)

% Ellipse parametrisation via eigendecomposition
[V, D] = eig(Cov2);
theta  = linspace(0, 2*pi, 200);
circle = [cos(theta); sin(theta)];
ellipse = V * sqrt(chiSq * D) * circle;

plot(mu(1) + ellipse(1,:), mu(2) + ellipse(2,:), ...
    '-', 'Color', color, 'LineWidth', 1.2)

end