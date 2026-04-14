function [M, P_pred1, P_pred2, PP] = compute_sensitivity_and_cov_derivatives(F, G, Sigma0, C, sigma2, u, K, ...
                                                                         F1, F2, G1, G2, Sigma01, Sigma02, Sigma, Sigma1, Sigma2)
% function [M, P_pred1, P_pred2, PP] = compute_sensitivity_and_cov_derivatives(F, G, Sigma0, C, sigma2, u, K, ...
%                                                                          F1, F2, G1, G2, Sigma01, Sigma02, Sigma, Sigma1, Sigma2)
% compute_sensitivity_and_cov_derivatives - Calcule récursivement
%   M_k = E[ dhat_{k|k-1}/dtheta1 * (dhat_{k|k-1}/dtheta2)' ]
%   ainsi que les dérivées des matrices de covariance de prédiction
%   P_pred1{k} = dP_{k|k-1}/dtheta1, P_pred2{k} = dP_{k|k-1}/dtheta2.
%
% Entrées :
%   F, G, Sigma0 : matrices du système discret (n x n, n x p, n x n)
%   C : matrice d'observation (q x n)
%   sigma2 : variance scalaire du bruit d'observation
%   u : matrice de commande (p x K)
%   K : horizon temporel
%   F1, F2 : dérivées de F par rapport à theta1, theta2 (n x n)
%   G1, G2 : dérivées de G (n x p)
%   Sigma01, Sigma02 : dérivées de Sigma0 (n x n)
%   Sigma : covariance stationnaire initiale (n x n)
%   Sigma1, Sigma2 : dérivées de Sigma (n x n)
%
% Sorties :
%   M : cellule de taille K, M{k} = matrice n x n (espérance du produit)
%   P_pred1 : cellule de taille K, P_pred1{k} = dP_{k|k-1}/dtheta1
%   P_pred2 : cellule de taille K, P_pred2{k} = dP_{k|k-1}/dtheta2
%   PP : cellule de taille K, P_{k|k-1}

    n = size(F,1);
    q = size(C,1);
    p = size(G,2);
    
    % Vérifier que u a la bonne dimension
    if size(u,2) < K
        error('u doit avoir au moins K colonnes');
    end
    if size(u,1) ~= p
        error('Le nombre de lignes de u doit correspondre à la deuxième dimension de G');
    end
    
    % Initialisation
    P_filt = Sigma;
    P_filt1 = Sigma1;
    P_filt2 = Sigma2;
    
    % Vecteur augmenté Z = [x̂; X1; X2] (3n x 1)
    EZ = zeros(3*n, 1);      % espérance (déterministe)
    CovZ = zeros(3*n);       % covariance (nulle initialement car Z0 déterministe)
    
    % Cellules de sortie
    M = cell(K,1);
    P_pred1 = cell(K,1);
    P_pred2 = cell(K,1);
    PP = cell(K,1);
    
    for k = 1:K
        % ---- Prédiction des covariances ----
        P_pred = F * P_filt * F' + Sigma0;
        PP{k} = P_pred;
        P_pred1{k} = F * P_filt1 * F' + F1 * P_filt * F' + F * P_filt * F1' + Sigma01;
        P_pred2{k} = F * P_filt2 * F' + F2 * P_filt * F' + F * P_filt * F2' + Sigma02;
        
        % ---- Innovation et gain ----
        S = C * P_pred * C' + sigma2 * eye(q);
        K_gain = (P_pred * C') / S;   % résout K_gain * S = P_pred * C'
        
        S1 = C * P_pred1{k} * C';
        S2 = C * P_pred2{k} * C';
        K1 = (P_pred1{k} * C' - K_gain * S1) / S;
        K2 = (P_pred2{k} * C' - K_gain * S2) / S;
        
        % ---- Mise à jour des covariances filtrées ----
        IKC = eye(n) - K_gain * C;
        P_filt_new = IKC * P_pred;
        P_filt1_new = IKC * P_pred1{k} - K1 * C * P_pred;
        P_filt2_new = IKC * P_pred2{k} - K2 * C * P_pred;
        
        % ---- Construction des matrices pour la mise à jour de Z ----
        % A (3n x 3n)
        A = zeros(3*n);
        A(1:n, 1:n) = F;
        A(n+1:2*n, 1:n) = IKC * F1;
        A(n+1:2*n, n+1:2*n) = IKC * F;
        A(2*n+1:3*n, 1:n) = IKC * F2;
        A(2*n+1:3*n, 2*n+1:3*n) = IKC * F;
        
        % B (3n x p)
        B = zeros(3*n, p);
        B(1:n, :) = G;
        B(n+1:2*n, :) = IKC * G1;
        B(2*n+1:3*n, :) = IKC * G2;
        
        % C_mat (3n x q)
        C_mat = zeros(3*n, q);
        C_mat(1:n, :) = K_gain;
        C_mat(n+1:2*n, :) = K1;
        C_mat(2*n+1:3*n, :) = K2;
        
        % ---- Mise à jour de la moyenne et covariance de Z ----
        EZ_new = A * EZ + B * u(:,k);
        CovZ_new = A * CovZ * A' + C_mat * S * C_mat';
        
        % ---- Calcul des sensibilités prédites Y = [X1_{k|k-1}; X2_{k|k-1}] ----
        H = zeros(2*n, 3*n);
        H(1:n, 1:n) = F1;
        H(1:n, n+1:2*n) = F;
        H(n+1:2*n, 1:n) = F2;
        H(n+1:2*n, 2*n+1:3*n) = F;
        
        D = [G1; G2];
        
        EY = H * EZ + D * u(:,k);
        CovY = H * CovZ * H';
        
        % ---- Matrice M_k ----
        % CovY est 2n x 2n, on prend le bloc (1:n, n+1:2*n)
        Cov12 = CovY(1:n, n+1:2*n);
        M{k} = Cov12 + EY(1:n) * EY(n+1:2*n)';
        
        % ---- Mise à jour pour le prochain pas ----
        EZ = EZ_new;
        CovZ = CovZ_new;
        P_filt = P_filt_new;
        P_filt1 = P_filt1_new;
        P_filt2 = P_filt2_new;
    end
end
