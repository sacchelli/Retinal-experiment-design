function [Y,Eta] = data_deterministic_1D(controls,idx,nn,sig2,...
                                    q, C, F, G, N)
% function [Y,Eta] = data_deterministic_1D(controls,idx,nn,sig2,...
%                                     q, C, F, G, N)
% idx = indices of inputs (elements of controls) that shoud be used, with
% nn_i>0 = nb. of repetitions of input = controls{idx(i)}
% Returns Y = mean of simulated data and Eta = mean responses for repeated
% observations at inputs defined by idx, both of size (q*N)*Lnn, 
% with Lnn = length(nn) the number of different inputs considered

% model response for input u -> eta(u)=H*u with
H = buildH(C, F, G, N);

Lnn=length(nn);
Ybar=zeros(q*N,Lnn);   % vectors of averaged observations
Eta=zeros(q*N,Lnn);    % vectors of model responses
for i_input=1:Lnn
    nb_repeat=nn(i_input);   
    % compute eta (deterministic response, length q*N) for i_input
    u = controls{idx(i_input)};
    u = u(:);
    Eta(:,i_input)=H*u;
    Y=zeros(q*N,nb_repeat);
    for i_repeat=1:nb_repeat
        % generate observations = eta + noise normal(0,sig2)
        Y(:,i_repeat)=Eta(:,i_input)+sqrt(sig2)*randn(N*q,1); 
    end
    Ybar(:,i_input)=mean(Y,2);       
end 
Y=Ybar;
end