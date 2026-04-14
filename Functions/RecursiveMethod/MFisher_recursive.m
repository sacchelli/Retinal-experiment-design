function [M] = MFisher_recursive(design_model,controls,theta, sig2,Delta, q, Acell, B, C,D, N)
% function [M] = MFisher_recursive(design_model,controls,theta, sig2,Delta, q, Acell, B, C,D, N)
% computes the Fisher information matrices for parameters theta and the
% input signals in controls by recursion
%                                ------------
% The calculation has a problem in the case design_model = 'stochastic' as
% the result does not coincide with that provided by MFisher.m which is
% based on the expansion of the response. As, morever, the calculations with
% MFisher_recursive.m are much slower than with MFisher.m for the case 
% design_model = 'stochastic', this function should NOT be used.
 
[K,~]=size(controls);
M = cell(K, 1);
p=length(theta);
[n,~]=size(Acell{1});

A = zeros(n,n);
for i = 1:p
    A = A+ theta(i) * Acell{i};
end

[F,G] = buildFG(A,B, Delta);
Sigma = lyap(A,D*D');           % Cov of initial state
Sigma0 = Sigma-F*Sigma*F';  % Cov of state noise

timebar= waitbar(0,'M Fisher computations...'); 
for i=1:K
    waitbar(i/K,timebar);
    if strcmp(design_model,'deterministic')
        R=zeros(q*N,p);
        for itheta=1:p
            [dEtadti] = der_response_deterministic_recursive1D(itheta,controls,i,1,...
                                    theta, Delta, prec, q, Acell, B, C, N);
            R(:,itheta)=dEtadti; % column vector (q*N)*1
        end
        M{i}=R'*R/(N*q)/sig2;
    elseif strcmp(design_model,'stochastic')
        % this is more difficult... we need a Kalman filter including derivatives etc.
        u = controls{i};
        Mi=zeros(p,p);
        Q=zeros(p,p);
        F_deriv=cell(p); 
        G_deriv=cell(p);
        Sigma_deriv=cell(p);
        Sigma0_deriv=cell(p);

        for j=1:p % (theta_j plays the role of theta_i)
                Ab = [A, zeros(n,n) 
                       Acell{j}, A];
                Bb = [B 
                       zeros(size(B))];
                [Fb,Gb] = buildFG(Ab,Bb, Delta);
                F_deriv{j}=Fb(n+1:2*n,1:n);
                G_deriv{j}=Gb(n+1:2*n,:);
                Sigmaj = lyap(A, Acell{j}*Sigma + Sigma*Acell{j}');
                Sigma_deriv{j}=Sigmaj;
                Sigma0_deriv{j} = Sigmaj - F_deriv{j} * Sigma * F' - F * Sigmaj * F' - F * Sigma * F_deriv{j}';
           
            for k=1:j
                % compute quantities depending on j and k with a Kalman filter
                [Exx, dPj, dPk, PP] = compute_sensitivity_and_cov_derivatives(F, G, Sigma0, C, sig2, u, N, ...
                                                                         F_deriv{j}, F_deriv{k}, G_deriv{j}, G_deriv{k}, ...
                                                                         Sigma0_deriv{j}, Sigma0_deriv{k}, Sigma, Sigma_deriv{j}, Sigma_deriv{k});
                % Element j,k of the FIM
                Mi_jk=0;
                for t=1:N                                       
                    djPt_tm1=dPj{t};
                    dkPt_tm1=dPk{t};
                    Pt_tm1=PP{t};

                    Omega_t=C*Pt_tm1*C'+sig2*eye(q); 
                    Omega_t_m1C=Omega_t\C;
                    Exxt=Exx{t};
                    %
                    Mi_jk=Mi_jk+ trAB(C'*Omega_t_m1C, Exxt) +...
                        0.5*trAB(Omega_t_m1C*djPt_tm1*C', Omega_t_m1C*dkPt_tm1*C');
                end
                Mi(j,k)=Mi_jk;
                Mi(k,j)=Mi_jk;
            end    
        end
        M{i}=Mi/(N*q);
    end
end    
close(timebar)

end