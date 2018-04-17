function [ R, alpha ] = fun_CLAML2( X,R0,lambda,mu )
%FUN_CLAML 此处显示有关此函数的摘要
%   此处显示详细说明
if nargin<3
    lambda = 3;
    mu = 1;
end
iter_num = 2;
[M,N] = size(X);
alpha = 0;
rou = 0;
R = eye(M);
for j = 1 : iter_num
    rou = 0;
    for i = 1 : N
        rou = rou + (X(:,i)' * inv(alpha * (R0-R) + R) * X(:,i) + 1/mu)\(X(:,i) * X(:,i)');
    end
    alpha = (lambda + M +1)/N * abs(trace((rou - R)/ (R0-R)))/M;
    alpha = max(min(1,alpha),0);
    R = alpha * (R0-R) + R;
end
end

