function [ R_AMLCC, alpha ] = fun_AMLCC(X,R_KA )
%%《Covariance matrixestimation for CFAR detection in correlated heavy tailed clutter》
%FUN_AML 此处显示有关此函数的摘要
%   此处显示详细说明
%近似最大似然协方差
%X:训练样本
[M,N] = size(X);
X_t = zeros(M,N);
R_AMLCC = eye(M,M);%R_x0;%eye(M,M);fun_NSCMN(X)%以单位阵为迭代初值是，第二次迭代结果为NSCM结果
alpha0 = 0;
% count = 1;
for i = 1:20%
    for j = 1:N
        tao = (1-alpha0) * (X(:,j)' * inv(R_AMLCC - alpha0*R_KA) * X(:,j))/M;
        X_t(:,j) =  X(:,j)/(tao^0.5);
    end
    [R_AMLCC,alpha] = fun_CC(X_t,R_AMLCC,R_KA );
    if abs(alpha0-alpha)<1e-2
        break;
    end
    alpha0 = alpha;
end
end

