function [Y,Eta,SigmaBoldHat] = data_stochastic_1D(controls,idx,nn,...
                                    q, C, F, G, SigmaBold, N)
% function [Y,Eta] = data_stochastic_1D(controls,idx,nn,sig2,...
%                                     q, C, F, G, N)
% idx = indices of inputs (elements of controls) that shoud be used, with
% nn_i>0 = nb. of repetitions of input = controls{idx(i)}
% Returns Y = mean of simulated data and Eta = mean responses for repeated
% observations at inputs defined by idx, both of size (q*N)*Lnn, 
% with Lnn = length(nn) the number of different inputs considered

% model response for input u -> eta(u)=H*u with
H = buildH(C, F, G, N);

Lnn=length(nn);
nnTot = sum(nn);
Ybar=zeros(q*N,Lnn);   % vectors of averaged observations
Eta=zeros(q*N,Lnn);    % vectors of model responses
SigmaBoldHat = zeros(q*N,q*N);
for i_input=1:Lnn
    nb_repeat=nn(i_input);
    % compute eta (deterministic response, length q*N) for i_input
    u = controls{idx(i_input)};
    u = u(:);
    Eta(:,i_input)=H*u;
    Yall = mvnrnd(Eta(:,i_input),SigmaBold,nb_repeat)';
    ybar = mean(Yall,2);
    SigmaBoldHat = SigmaBoldHat + nb_repeat/nnTot*(Yall*Yall'/nb_repeat - ybar*ybar');
    Ybar(:,i_input)=ybar;     
end 
Y=Ybar;
end