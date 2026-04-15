function Jtrace = trace_part_JML(SigmaBoldInv, nn, CovY)
% function Jtrace = trace_part_JML(SigmaBoldInv, nn, CovY)
% compute the trace part (due to the use of a sufficient statistic) of the ML criterion for a stochastic model

Nexp=sum(nn);
BigCov=zeros(size(CovY{1}));
for i=1:length(nn)
    BigCov=BigCov+(nn(i)/Nexp)*CovY{i};
end
Jtrace=trAB(SigmaBoldInv,BigCov);
end