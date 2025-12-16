addpath './Functions'
clc

tic


rng(1)

%%%% Data // following the reference.

delta = 0.00389 % cm
sigmaBG = 0.15 % cm
sigmaAG = 0.15 % cm



%%%% Numerical data

dn = 10; % size of one layer

n = 3 * dn;

m = dn;

p = 7;

q = dn;

N = 100; % Time window


%%%% Model

Gamma = laplacian1DGraph(dn);
GammaBG = gaussianPooling1DGraph(dn,delta,sigmaBG);

