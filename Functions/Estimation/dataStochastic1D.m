function [Y, Eta, SigmaBoldHat] = dataStochastic1D(controls, idx, nn, ...
                                    q, C, F, G, SigmaBold, N)
% function [Y, Eta, SigmaBoldHat] = dataStochastic1D(controls, idx, nn,
%                                     q, C, F, G, SigmaBold, N)
% idx = indices of inputs (elements of controls) that should be used, with
% nn_i > 0 = nb. of repetitions of input = controls{idx(i)}
% Returns Y = mean of simulated data and Eta = mean responses for repeated
% observations at inputs defined by idx, both of size (q*N)*Lnn,
% with Lnn = length(nn) the number of different inputs considered.
% Model response for input u -> eta(u) = H*u with

H = buildH(C, F, G, N);

Lnn = length(nn);
nnTot = sum(nn);

Ybar = zeros(q*N, Lnn);        % vectors of averaged observations
Eta = zeros(q*N, Lnn);         % vectors of model responses
SigmaBoldHat = zeros(q*N, q*N);

for iInput = 1:Lnn

    nRepeat = nn(iInput);

    % Compute eta (deterministic response, length q*N) for iInput.
    u = controls{idx(iInput)};
    u = u(:);
    Eta(:, iInput) = H * u;

    Yall = mvnrnd(Eta(:, iInput), SigmaBold, nRepeat)';
    yBar = mean(Yall, 2);
    SigmaBoldHat = SigmaBoldHat + nRepeat/nnTot * (Yall*Yall'/nRepeat - yBar*yBar');
    Ybar(:, iInput) = yBar;

end

Y = Ybar;

end