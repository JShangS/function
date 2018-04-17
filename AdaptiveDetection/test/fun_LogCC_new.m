function [ R,alpha ] = fun_LogCC_new( X,R_KA );
%FUN_LOGCC_NEW 此处显示有关此函数的摘要
%   此处显示详细说明
[N,L] = size(X);
Ri = zeros(N,N,L);
logm_R = zeros(N,N);
for i = 1:L
    Ri(:,:,i) = fun_Positive(X(:,i),3);
    logm_R = logm_R + fun_Logm(Ri(:,:,i))/L;
end
t1 = 0;
t2 = 0;
for i = 1:L
    t11 = fun_Logm(Ri(:,:,i)) - logm_R;
    t1 = t1 + norm(t11,'fro')^2;
    t22 = fun_Logm(R_KA) - fun_Logm(Ri(:,:,i));
    t2 = t2 + norm(t22,'fro')^2;
end
alpha = t1/t2;
R = alpha * fun_Logm(R_KA) + (1-alpha)*logm_R;
R = fun_Expm(R);