function [ R_AML, alpha ] = fun_AMLCC2( X,R_KA )
%%《Covariance matrixestimation for CFAR detection in correlated heavy tailed clutter》
%FUN_AML 此处显示有关此函数的摘要
%   此处显示详细说明
%近似最大似然协方差
%X:训练样本
[M,N] = size(X);
R_AML = eye(M,M);%R_x0;%eye(N,N);%以单位阵为迭代初值是，第二次迭代结果为NSCM结果
R_AML0 = zeros(M,M);
% count = 1;
alpha0 = 1;
alpha = 0;
for j = 1:10
%     R_AML_inv = inv(R_AML);
    alpha0=alpha;
    R_AML0 = R_AML;
%     tao_child = diag((X'/R_AML0*X)/M)；
    t1 = 0;
    t2 = 0;
    t3 = 0;
    for i = 1:N
        tao = X(:,i)'/R_AML0*X(:,i)/M;
        t1 = t1 + norm(X(:,i) * X(:,i)'/tao-R_AML0,'fro')^2;
        t2 = t2 + norm(R_KA - X(:,i) * X(:,i)'/tao,'fro')^2;
        t3 = t3 + X(:,i) * X(:,i)'/N/tao;
    end
    alpha = max(min(1,t1/t2),0);
    R_AML = alpha * R_KA + (1-alpha) * t3;
    if abs(alpha - alpha0)<1e-2
        break;
    end
%     R_AML = t3;
end
end

