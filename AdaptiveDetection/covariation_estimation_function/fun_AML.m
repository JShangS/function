function [ R_AML ] = fun_AML( X )
%%《Covariance matrixestimation for CFAR detection in correlated heavy tailed clutter》
%FUN_AML 此处显示有关此函数的摘要
%   此处显示详细说明
%近似最大似然协方差
%X:训练样本
[M,N] = size(X);
R_AML = eye(M,M);%R_x0;%eye(N,N);%以单位阵为迭代初值是，第二次迭代结果为NSCM结果
tao_child = 1;%%本次迭代值
tao_parent = 0;%%上次迭代值
% count = 1;
while (abs(tao_child-tao_parent)>0.1)%
    tao_parent = tao_child;
    R_AML_inv = inv(R_AML);
    tao_child = abs(X'*R_AML_inv*X)/M;
    R_AML_t = 0;
    for i = 1:N
        R_AML_t = R_AML_t+X(:,i)*X(:,i)'/N;
    end
    R_AML = abs(R_AML_t);
%     count =count+1;
%     if count >=3  
%         break;
%     end
end
end

