function C = blkadd(A, B)

C = [A, zeros(size(A,1),size(B,2)); zeros(size(B,1),size(A,2)), B];
end