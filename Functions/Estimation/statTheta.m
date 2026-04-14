function [meanTheta, meanAbsBiasTheta, covTheta, sqrtEigCov] = statTheta(THETA, thetaTrue)
% function [meanTheta, meanAbsBiasTheta, covTheta, sqrtEigCov] = statTheta(THETA, thetaTrue)
%
% THETA is a N*p array: p parameters, N row vectors of estimated values theta_i. Returns:
% meanTheta = mean(theta_i),
% meanAbsBiasTheta = mean(|theta_i - thetaTrue|),
% covTheta = covariance matrix (p*p) of the theta_i,
% sqrtEigCov = sqrt of the eigenvalues of covTheta.

[N, p] = size(THETA);

covTheta = cov(THETA);
meanTheta = mean(THETA, 1);
meanAbsBiasTheta = mean(abs(THETA - ones(N,1)*thetaTrue), 1);
sqrtEigCov = sqrt(eig(covTheta)');

end