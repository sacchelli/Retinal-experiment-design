function [Cov_theta_theoretical,sqrt_eig_Cov_theoretical] = Stat_THETA_theoretical(M,idx,nn,N,q)
% 
% 
% Cov_theta_theoretical = theoretical (i.e., asymptotic) covariance matrix (p*p) 
%   of the ML estimator (= inverse of the Fisher information matrix)
% sqrt_eig_Cov_theoretical = sqrt of its eigenvalues

Lnn=length(nn);

MFisher=0*M{1};
for i_input=1:Lnn
    nb_repeat=nn(i_input);   
    MFisher = MFisher + nb_repeat*M{idx(i_input)};
end    
Cov_theta_theoretical=inv(MFisher)/(N*q);
sqrt_eig_Cov_theoretical = sqrt(eig(Cov_theta_theoretical)');
end