function value = trMinv(M)
    if rank(M) < size(M, 1)
        value = Inf;
    else
        value = trace(inv(M));
    end
end