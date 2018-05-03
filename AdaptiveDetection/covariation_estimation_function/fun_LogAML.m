function [ R_logAML ] = fun_LogAML( X )
%FUN_LOGAML 此处显示有关此函数的摘要
%   此处显示详细说明
[M,N] = size(X);
R_logAML = eye(M,M);%R_x0;%eye(N,N);%以单位阵为迭代初值是，第二次迭代结果为NSCM结果
tao_child = 1;%%本次迭代值
tao_parent = 0;%%上次迭代值
% count = 1;
while (abs(tao_child-tao_parent)>0.5)%
    tao_parent = tao_child;
    iR_AML = inv(R_logAML);
    tao_child = diag((X'*iR_AML*X)/M);
    R_AML_t = 0;
    for i = 1:N
%         tao = (X(:,i)'*R_AML_inv*X(:,i))/M;
        Ri = fun_Positive(X(:,i),4);
        R_AML_t = R_AML_t+logm(Ri/tao_child(i));
    end
    R_AML_t = expm(R_AML_t/N);
%     if norm(R_AML_t-R_logAML,'fro')<1e-2
%         break;
%     end
    R_logAML = R_AML_t;
end

end

