function H = buildH(C, F, G, N)

n = size(C,2);
q = size(C,1);
m = size(G,2);

% dHk = eye(size(F));
% 
% col = [C*dHk*G];
% 
% for i=2:N
%     dHk=dHk*F;
%     col=[col;C*dHk*G];
% end
% 
% mat=col;
% 
% for i=2:N
%     col=[zeros(q,m);col(1:end-q,:)];
%     mat=[mat,col];
% end


dHk = eye(size(F));
col = zeros(N*q,m);
col(1:q,:) = C*dHk*G;
for i=2:N
    dHk=dHk*F;
    col((i-1)*q+1:i*q,:)=C*dHk*G;
end

mat = zeros(N*q,N*m);

for i=1:N
    mat((i-1)*q+1:end,(i-1)*q+1:i*q) = col(1:(N-i+1)*q,:);
end

H=mat;

end