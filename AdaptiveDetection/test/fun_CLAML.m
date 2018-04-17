function [ R, alpha ] = fun_CLAML( X,R0,R1,lambda,mu )
%FUN_CLAML 此处显示有关此函数的摘要
%   此处显示详细说明
if nargin<4
    lambda = 3;
    mu = 1;
end
iter_num = 3;
[M,N] = size(X);
alpha = 0.01;
% R01 = fun_Positive(R0-R1);
R01=R0-R1;
for j = 1 : iter_num
    rou = 0;
    for i = 1 : N
        rou = rou + ( trace((alpha * (R01) + R1)\(X(:,i) * X(:,i)')) + 1/mu)\(X(:,i) * X(:,i)');
    end
    alpha = (lambda + M +1)/N * abs(trace((rou -R1)/ (R01)))/M;
    alpha = max(min(1,alpha),0);
end
R = alpha * (R0-R1) + R1;
end

