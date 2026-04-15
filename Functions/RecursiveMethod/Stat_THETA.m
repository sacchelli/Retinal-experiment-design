function [mean_theta,mean_absbias_theta,Cov_theta,sqrt_eig_Cov] = Stat_THETA(THETA,theta_true)
% function [mean_theta,mean_absbias_theta,Cov_theta,sqrt_eig_Cov] = Stat_THETA(THETA,theta_true)
% 
% THETA is a N*p array : p parameters, N row vectors of estimated values theta_i; returns:
% mean_theta = mean(theta_i), 
% mean_absbias_theta = mean(|theta_i-theta_true|),
% Cov_theta = covariance matrix (p*p) of the theta_i
% sqrt_eig_Cov = sqrt of the eigenvalues of Cov_theta

[N,p]=size(THETA);

Cov_theta=cov(THETA);
mean_theta=mean(THETA,1);

mean_absbias_theta=mean( abs(THETA-ones(N,1)*theta_true) ,1);

sqrt_eig_Cov=sqrt(eig(Cov_theta)');

end