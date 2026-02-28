function [nn,indices] = efficient_rounding(w,n)
% function [nn,indices] = efficient_rounding(w,n)
%   w a vector of L weights w_i >= 0 
%   n an integer (can be larger or smaller than L)
% The function returns
% a row vector nn integers nn_i such that nn_i/n approximates w_i
% a row vector of indices indicating which i correspond to the w_i used
% Based on the efficient rounding method of [Pukelsheim & Reider, Biometrika, 1992].
% If L<= n, an nn_i>0 is associated with each w_i (i.e., length(nn)=L)
% otherwise, small weights are neglected and length(nn) may be less than n 

w=w(:)';
w=w/sum(w);
L=length(w);

if L<n
    ibig=1:L;
    wbig=w;
    Lbig=L;
else
    % We want to remove small weights that will be approximated by zero in any case
    [wsorted, isorted]=sort(w,'descend'); % sort the weights by decreasing values
    % --> find a threshold tau such that only w(i)>tau are considered
    %       tau = wsorted(icut) for some icut
    %   a) for sure, only the n largest weights should be used
        icut=n;
    %   b) we can cut more than that: 
    % %       any wsorted(i+1)< wsorted(i)/nn_i is useless, and necessarily nn_i <= ceil(n/i)
    %     wsorted_n=wsorted(1:icut); 
    %     sum_wsorted_n=sum(wsorted_n); 
    %     Ln=length(wsorted_n);
    %     wsorted_n=wsorted_n/sum_wsorted_n;
    %     %wsorted(1:Ln)
    %     sequence=wsorted_n-[0 wsorted_n(1:Ln-1)./ceil(n.*(1:Ln-1).^(-1))];
    %     % (this is more efficient than sequence=wsorted_n-[0 wsorted_n(1:Ln-1)]/n;)
    %     ii=find(sequence<0);
    %     ii_1=n;
    %     if ~isnan(ii)
    %         ii_1=ii(1);
    %     end
    %     icut=min(icut,ii_1);
    % %   now cut 
    %     threshold=wsorted_n(icut);

        wsorted_n=wsorted(1:icut);
        m=ceil(n*wsorted_n); % --> this yields sum(m)>n points
        while sum(m)>n
            % remove some
            [~,imax]=max(m/n-wsorted_n);
            m(imax)=m(imax)-1;
        end

        % remove points with m(i)=0
        iremove=find(m==0);
        if isempty(iremove)==0
            icut=iremove(1)-1; % --> smallest weight to be kept
        end        
        % threshold=wsorted_n(icut);
        % ibig=find(w>=threshold);
        ibig=isorted(1:icut);
        wbig=w(ibig);
        Lbig=length(wbig);
end  
% Now, the efficient rounding method of [Pukelsheim & Reider, Biometrika, 1992]
%   first step: allocate initial numbers to big enough weigths    
nnbig=ceil((n-Lbig/2)*wbig);

%   second step: adjust a few allocations
while sum(nnbig)<n
    [~,istar]=min(nnbig./wbig);
    nnbig(istar)=nnbig(istar)+1;
end
while sum(nnbig)>n
    [~,istar]=max((nnbig-1)./wbig);
    nnbig(istar)=nnbig(istar)-1;
end
indices=ibig;
nn=nnbig;

end

