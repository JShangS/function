function [ distance ] = fun_KLD( A,B )
%%Geometric means and medians with applications to target detection(12)ʽ
[N,~] = size(A);
distance = trace(inv(B) * A - eye(N)) - log(det(inv(B) * A));
end