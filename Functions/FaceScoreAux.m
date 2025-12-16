function [score, weights] = FaceScoreAux(M)

% Je fais l'hypothèse que M est de la forme p x p x K

p = size(M,1);
K = size(M,3);
         
step = 0.01;        % resolution
vals = 0:step:1;   % possible weight values

combs = [];        % store results

for w1 = vals
    for w2 = vals
        w3 = 1 - w1 - w2; % enforce sum = 1
        if w3 >= 0 && w3 <= 1
            combs = [combs; [w1, w2, w3]];
        end
    end
end

maxScore = -Inf;
wStar=combs(1,:);

for i = 1:size(combs,1)
    w = combs(i,:);
    Mmu = sum(M .* reshape(w, 1, 1, K), 3);
    MmuScore = log(det(Mmu));
    if MmuScore > maxScore
        wStar = w;
        maxScore = MmuScore;
    end
end


score = maxScore;
weights = wStar;


end