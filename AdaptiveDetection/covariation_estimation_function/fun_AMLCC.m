function [ R_AML, alpha ] = fun_AMLCC(X,R_KA )
%%《Covariance matrixestimation for CFAR detection in correlated heavy tailed clutter》
%FUN_AML 此处显示有关此函数的摘要
%   此处显示详细说明
%近似最大似然协方差
%X:训练样本
[M,N] = size(X);
X_t = zeros(M,N);
R_AML = fun_AML(X);%R_x0;%eye(M,M);fun_NSCMN(X)%以单位阵为迭代初值是，第二次迭代结果为NSCM结果
[R_AML,alpha0] = fun_CC(X,R_AML,R_KA );
alpha = 0;
% count = 1;
for i = 1:20%
%     i
%     tao_child = diag((X'/R_AML*X)/M);
    iR_AML = inv(R_AML);
    for j = 1:N
        tao = (X(:,j)'*X(:,j))/M;
        X_t(:,j) = X(:,j)/(tao)^(1/2);
    end
    [R_AML,alpha] = fun_CC(X_t,R_AML,R_KA );
    if abs(alpha0-alpha)<1e-2
        break;
    end
    alpha0 = alpha;
end
end

