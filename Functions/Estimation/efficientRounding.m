function [nn, indices] = efficientRounding(w, n)
% function [nn, indices] = efficientRounding(w, n)
%   w a vector of L weights w_i >= 0
%   n an integer (can be larger or smaller than L)
% The function returns
% a row vector nn of integers nn_i such that nn_i/n approximates w_i
% a row vector of indices indicating which i correspond to the w_i used
% Based on the efficient rounding method of [Pukelsheim & Reider, Biometrika, 1992].
% If L <= n, an nn_i > 0 is associated with each w_i (i.e., length(nn) = L)
% otherwise, small weights are neglected and length(nn) may be less than n

w = w(:)';
w = w/sum(w);
L = length(w);

if L < n

    iBig = 1:L;
    wBig = w;
    LBig = L;

else

    % We want to remove small weights that will be approximated by zero in any case.

    [wSorted, iSorted] = sort(w, 'descend'); % sort the weights by decreasing values
    
    iCut = n;

    wSortedN = wSorted(1:iCut);
    m = ceil(n * wSortedN); % --> this yields sum(m) > n points

    while sum(m) > n
        % remove some
        [~, iMax] = max(m/n - wSortedN);
        m(iMax) = m(iMax) - 1;
    end

    % remove points with m(i) = 0
    iRemove = find(m == 0);

    if isempty(iRemove) == 0
        iCut = iRemove(1) - 1; % --> smallest weight to be kept
    end

    iBig = iSorted(1:iCut);
    wBig = w(iBig);
    LBig = length(wBig);

end

% Now, the efficient rounding method of [Pukelsheim & Reider, Biometrika, 1992].

%   First step: allocate initial numbers to big enough weights.
nnBig = ceil((n - LBig/2) * wBig);

%   Second step: adjust a few allocations.
while sum(nnBig) < n
    [~, iStar] = min(nnBig ./ wBig);
    nnBig(iStar) = nnBig(iStar) + 1;
end

while sum(nnBig) > n
    [~, iStar] = max((nnBig - 1) ./ wBig);
    nnBig(iStar) = nnBig(iStar) - 1;
end

indices = iBig;
nn = nnBig;

end