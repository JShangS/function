function [ R,w] = fun_DWCM( X,R_KA,x0 )
%FUN_DWCM 此处显示有关此函数的摘要
%   此处显示详细说明
%%先随便取个名字，距离度量的权值，协方差估计
%%X:训练数据
%%RKA:先验协方差
%%x0:检测单元
[M,N] = size(X);
distance = zeros(N,1);  %%距离度量
Sx0 = x0*x0';
for i = 1:N
    distance(i) = norm(X(:,i)*X(:,i)'-Sx0,'fro');
end
distance(N+1) = norm(R_KA - Sx0,'fro');
w = exp(-distance)./sum(exp(-distance));
R_t = 0;
for i =1:N
    R_t = R_t + w(i)*(X(:,i)*X(:,i)');
end
R_t = R_t + w(N+1)*R_KA;
R_Child = R_t;
R_Parent = zeros(size(R_Child));
count = 1;
while(sum(sum(abs(R_Parent-R_Child)))>0.001)
    sprintf('%d',count);
    R_Parent = R_Child;
    for i = 1:N
        distance(i) = norm(X(:,i)*X(:,i)'-R_Child,'fro');
    end
    distance(N+1) = norm(R_KA - R_Child,'fro');
    w = exp(-distance)./sum(exp(-distance));
    R_t = 0;
    for i =1:N
        R_t = R_t + w(i)*(X(:,i)*X(:,i)');
    end 
    R_t = R_t + w(N+1)*R_KA;
    R_Child = R_t;
    count = count+1;
end
R = R_Child;
end

