function [ R_AML, alpha0 ] = fun_AMLCC3(X,R_KA )
%%《Covariance matrixestimation for CFAR detection in correlated heavy tailed clutter》
%FUN_AML 此处显示有关此函数的摘要
%   此处显示详细说明
%近似最大似然协方差
%X:训练样本
[M,N] = size(X);
 R_AML = fun_NSCMN(X);
[R_AML,alpha0] = fun_CC(X,R_AML,R_KA );
for i = 1:20%
    [R_AML,alpha] = fun_CC(X,R_AML,R_KA );
    if abs(alpha0-alpha)<1e-2
        break;
    end
    alpha0 = alpha;
end
end

