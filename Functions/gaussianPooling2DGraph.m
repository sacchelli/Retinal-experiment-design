function M = gaussianPooling2DGraph(dn,delta,sigma)
    % 
    %  dn: number of points per dimension, one layer needs to be dn x dn.
    %   
    
    e1 = exp(-delta^2/(2*sigma^2));
    M = zeros(dn^2,dn^2);
    for k = 1:dn^2
        ik = mod(k,dn);
        jk = floor((k-1)/dn)+1;
        for l = 1:dn^2
            il = mod(l,dn);
            jl = floor((l-1)/dn)+1;
            M(l,k) = e1^((ik-il)^2+(jk-jl)^2);
        end
    end
end