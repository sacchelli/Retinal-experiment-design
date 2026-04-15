function logL = kalman_loglikelihood(Y, u, F, B, Sigma0, C, Sigmap, Sigma)
% KALMAN_LOGLIKELIHOOD Calcule la log-vraisemblance d'un modèle d'espace d'état linéaire gaussien
%   logL = kalman_loglikelihood(Y, u, F, B, Sigma0, C, Sigmap, Sigma)
%
%   Y  : matrice des observations (dimension q x N)
%   u  : matrice des commandes (dimension p x N). Si pas de commande, mettre zeros(p,N) ou []
%   F  : matrice de transition (n x n)
%   B  : matrice d'entrée (n x p)
%   Sigma0 : covariance du bruit d'état (n x n)
%   C  : matrice d'observation (q x n)
%   Sigmap : covariance du bruit d'observation (q x q)
%   Sigma : covariance initiale de l'état (n x n) (généralement la solution de l'équation de Lyapunov)
%
%   La log-vraisemblance est calculée par le filtre de Kalman.

    % Vérifications des dimensions
    [q, N] = size(Y);
    n = size(F,1);
    if size(F,2) ~= n
        error('F doit être carrée');
    end
    if size(B,1) ~= n
        error('B doit avoir n lignes');
    end
    p = size(B,2);
    if ~isempty(u)
        if size(u,2) ~= N
            error('u doit avoir N colonnes');
        end
        if size(u,1) ~= p
            error('Le nombre de lignes de u doit correspondre au nombre de colonnes de B');
        end
    else
        u = zeros(p, N); % si u non fourni, on met des zéros
    end
    
    if size(C,1) ~= q || size(C,2) ~= n
        error('C doit être de dimension q x n');
    end
    if any(size(Sigmap) ~= q)
        error('Sigmap doit être q x q');
    end
    if any(size(Sigma0) ~= n)
        error('Sigma0 doit être n x n');
    end
    if any(size(Sigma) ~= n)
        error('Sigma doit être n x n');
    end
    
    % Initialisation
    x_pred = zeros(n,1);
    P_pred = Sigma;
    
    logL = 0;
    for k = 1:N
        % Innovation
        y_pred = C * x_pred;
        innov = Y(:,k) - y_pred;
        S = C * P_pred * C' + Sigmap;
        
        % Calcul de la log-vraisemblance (décomposition de Cholesky pour stabilité)
        [L, p] = chol(S, 'lower');
        if p ~= 0
            warning('Matrice de covariance S non définie positive, ajout d''une régularisation');
            S = S + 1e-6 * eye(q);
            [L, p] = chol(S, 'lower');
            if p ~= 0
                error('S toujours non définie positive');
            end
        end
        logdetS = 2 * sum(log(diag(L)));
        innov_norm = L \ innov; % résout L * x = innov
        quad = innov_norm' * innov_norm;
        
        logL = logL - 0.5 * (logdetS + quad);
        
        % Gain de Kalman
        K = (P_pred * C') / S;  % résout K*S = P_pred*C'
        
        % Mise à jour
        x_upd = x_pred + K * innov;
        P_upd = (eye(n) - K * C) * P_pred;
        
        % Prédiction pour le prochain pas
        if k < N
            x_pred = F * x_upd + B * u(:, k+1);
            P_pred = F * P_upd * F' + Sigma0;
        end
    end
    
    % Constante
    %logL = logL - (N * q / 2) * log(2 * pi);
end
