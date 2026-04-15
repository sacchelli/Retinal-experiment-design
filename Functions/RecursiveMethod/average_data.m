function [Ymean,CovY] = average_data(Y,nn)
% function [Ymean,CovY] = average_data(Y,nn)
% Computes a sufficient statistic fof the data Y that contains repetitions
%   under the same experimental conditions
% Ymean = mean of simulated data for repeated observations 
% Y(:,i) = (N*q)*1 vector of observations at inputs nb. i, repeated nn(i) times,
% Ymean=zeros(q*N,Lnn) with Lnn = length(nn) the number of different inputs
% Ymean(:,i) = (N*q)*1 vector of everaged observations at input nb. i
% CovY = cell(Lnn, 1); each cell has size (q*N)*(q*N) and contains
%   the covariance matrix of the Y(:,i) for the nn(i) repetitions with input nb. i,
%   i=1,...,Lnn (the cell nb. i equals zeros(q*N,q*N) if no repetition for
%   this i)

[qN,Nexp]=size(Y);
if sum(nn)~=Nexp           % total number of experiments
    disp('Y should have sum(nn) columns')
    Ymean=[]; CovY=[];
    return
end    
Lnn=length(nn);
CovY = cell(Nexp, 1);
Ymean=zeros(qN,Lnn);  % vectors of averaged observations
first_index_Y=1;
for i_input=1:Lnn
    nb_repeat=nn(i_input);  
    Y_i_input=Y(:,first_index_Y:first_index_Y+nb_repeat-1);
    Ymean(:,i_input)=mean(Y_i_input,2);
    CovY{i_input}=zeros(qN,qN);
    if nb_repeat>1
        CovY{i_input}=cov(Y_i_input',1);
    end
    first_index_Y=first_index_Y+nb_repeat;
end 
end